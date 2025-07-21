function [lawn] = analyzeLawnProfiles(image,background,template,lawnRadius,scale)
% [lawn] = ANALYZELAWNPROFILES(image,background,template,lawnRadius,scale)
%
%   ANALYZELAWNPROFILES takes in an image of fluorescent bacterial patches
%   and normalizes that image to its background to account for lighting
%   inconsistencies. Lawns are detected using GETLAWNSTHRESHOLD and
%   intenstiy profiles are extracted and analyzed.
%
%   INPUTS:
%       - image [MxN num]: fluorescent image of bacterial patches
%           expressing GFP
%       - background [MxN double]: image of the filtered/fit light source
%       - template [MxN logical]: mask of the area to analyze
%       - lawnRadius [double]: expected radius of the lawns in pixels
%       - scale [double]: scale of the images in pixels/mm
%
%   OUTPUTS:
%       - lawn [struct]: structure of analysis
%           - mask [MxN logical]: mask of lawns dilated to ~0.5 mm from
%               the lawn edges
%           - labeled [MxN double]: labeled mask of dilated lawns
%           - image [MxN num]: original image
%           - centers [numLawnsx2 double]: X and Y positions of the center
%               of each labeled lawn
%           - radii [numLawnsx1 double]: radius of each labeled lawn
%               (before dilation) in pixels
%           - circularity [numLawnsx1 double]: circularity of each labeled
%               lawn (before dilation)
%           - imageNormalized [MxN num]: image divided by the background
%           - distances [1xnumBins double]: distances from the edge of each
%               lawn in mm; step size is 1e-3 (i.e. 1 um); negative values
%               are outside the lawn, positive values are inside
%           - pixelValuesNormalized [nmLawnsxnumBins double]: interpolated
%               pixel intensity values binned at each distance from the
%               patch edge
%           - edge [MxN logical]: mask of lawn edge
%           - analysis [numLawnsx14 table]: table of analysis variables
%               - xPeak: distance from edge to max pixel intensity (mm)
%               - yPeak: max pixel intensity
%               - xOuterEdge: distance of patch edge (0 mm)
%               - yOuterEdge: pixel intensity at patch edge
%               - xHalfMaxOuter: distance at which the pixel intensity has
%                   risen half way to maximum (mm)
%               - xHalfMaxInner: distance at which the pixel intensity has
%                   fallen half way to maximum (mm)
%               - yHalfMax: pixel intensity at half maximum
%               - FWHM: xHalfMaxInner - xHalfMaxOuter (mm)
%               - lawnRadius: radius of the lawn as computed from the area
%                   of the edge mask (mm)
%               - circularity: circularity of each lawn
%               - borderAmplitude: yPeak - yOuterEdge
%               - meanAmplitude: average pixel intensity as computed from
%                   the edge mask applied to imageNormalized
%               - centerAmplitude: average pixel intensity of the middle of
%                   the lawn as computed from the edge mask dilated by 50%
%                   applied to imageNormalized
%               - borderCenterRatio: borderAmplitude/centerAmplitude
%    
%   Written 2/21/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEGFP, GETLAWNSTHRESHOLD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Replace outliers (saturated pixels)
outlier = image == intmax(class(image));
imageOutlier = double(image);
imageOutlier(outlier) = NaN;
imageOutlier = fillmissing(imageOutlier,'linear');

% Normalize image to light source
imageNormalized = max(background,[],'all').*imageOutlier./background;

% Invert and blur image
imageInverse = imcomplement(imageNormalized);

% Detect patches via thresholding
[lawn] = getLawnsThreshold(imageInverse,template,lawnRadius);
lawnObjects = regionprops(lawn.labeled,'Circularity');
lawn.circularity = [lawnObjects.Circularity]';

% Dilate patch masks by 0.5 mm
lawn.image = image;
lawn.imageNormalized = imageNormalized;
lawn.mask = bwdist(lawn.mask) <= 0.5*scale;
lawn.labeled = bwlabel(lawn.mask);
lawnObjects = regionprops(lawn.labeled,'Area');
numLawn = length(unique(lawn.labeled(:)))-1;

% If fluorescence signal was not strong enough to find patches, skip
if numLawn == 0
    lawn.distances = [];
    lawn.pixelValuesNormalized = [];
    lawn.edge = lawn.mask;
    lawn.analysis = table();
else
    % Create distance bins that will be spaced such that each contains approx.
    % the same number of pixels
    binArea = 1e-2; % each bin will average this area (mm^2) worth of pixels
    r0 = sqrt(binArea/pi); % smallest bin radius
    rMax = ceil((sqrt(max([lawnObjects.Area])./pi))/scale + 0.1); % largest bin radius
    numBins = ceil((rMax/r0)^2);
    distances = r0*sqrt(0:numBins);

    % Upsample images via nearest neighbor interpolation
    maskUp = imresize(bwdist(~lawn.mask),3);
    labeledUp = bwlabel(round(imresize(lawn.mask,3)));
    imageUp = imresize(lawn.image,3);
    imageNormUp = imresize(lawn.imageNormalized,3);

    % Get average pixel values in each distance bin
    pixelValues = nan(numLawn,length(distances));
    pixelValuesNormalized = nan(numLawn,length(distances));
    for numROI = 1:numLawn
        maskROI = maskUp(labeledUp == numROI);
        [~,~,loc] = histcounts(-maskROI + max(maskROI),distances.*scale);
        pixelValues(numROI,unique(loc)) = ...
            groupsummary(imageUp(labeledUp == numROI),loc,'mean');
        pixelValuesNormalized(numROI,unique(loc)) = ...
            groupsummary(imageNormUp(labeledUp == numROI),loc,'mean');
    end

    % Interpolate distances and pixel values for evenly spaced bins (1 um)
    xMax = max(distances(sum(~isnan(pixelValuesNormalized),1) > 0));
    xStep = 1e-3;
    xInterp = 0:xStep:ceil(xMax/xStep)*xStep;
    yInterp = nan(numLawn,size(xInterp,2));
    for numROI = 1:numLawn
        yInterp(numROI,:) = interp1(distances,pixelValuesNormalized(numROI,:),xInterp);
    end

    xLag = zeros(size(yInterp,1),1);
    for j = 1:numLawn
        ySpline = csaps(xInterp,rescale(yInterp(j,:)),1-5e-6,xInterp);

        % Calculate curvature
        dy = gradient(ySpline,xInterp);
        dy2 = gradient(dy,xInterp);
        kappa = dy2./sqrt((1+dy.^2).^3);
        lawn.kappa(j,:) = kappa;

        % Use curvature to find outer patch edge
        [kappaVals,kappaPeaks,~,kappaProms] = findpeaks(log(max(kappa,0)),xInterp);
        dy(ySpline < 0.2) = NaN;
        [~,indRise] = min(dy);
        maxVal = max(kappaVals(kappaPeaks > xInterp(indRise) & isinf(kappaProms)));
        if xInterp(indRise) > 2
            xLag(j) = kappaPeaks(kappaVals == maxVal);
        else
            xLag(j) = kappaPeaks(find(kappaPeaks > xInterp(indRise) & ...
                isinf(kappaProms),1,'first')); % kappaPeaks > prctile(kappaPeaks,70)
        end
    end
    xOuterEdge = min(xLag);
    indLag = round((xLag - xOuterEdge)/xStep);

    % Shift arrays by offset
    xInterp = [-max(indLag)*xStep:xStep:-xStep,xInterp];
    yInterp = padarray(yInterp,[0 max(indLag)],NaN,'pre');
    yShift = nan(size(yInterp));
    for numROI = 1:size(yInterp,1)
        yShift(numROI,:) = circshift(yInterp(numROI,:),-indLag(numROI));
    end

    % Save lawn distances and pixel values relative to patch edge
    lawn.distances = flip(-(xInterp - xOuterEdge));
    lawn.pixelValuesNormalized = flip(yShift,2);

    % For each patch, find the inner and outer edge of the patch border as well
    % as the border peak and use these to compute metrics
    lawn.edge = false(size(lawn.mask));
    lawn.analysis = table();
    for j = 1:numLawn

        % Get distances and pixel values for this lawn (excluding NaN)
        indFit = ~isnan(lawn.pixelValuesNormalized(j,:));
        x = lawn.distances(indFit);
        y = lawn.pixelValuesNormalized(j,indFit);

        % Get peak of pixel values
        xPeak = interp1(y,x,max(y(x<1)));
        yPeak = interp1(x,y,xPeak);

        % Get pixel intensity at outer edge
        xOuterEdge = 0;
        yOuterEdge = interp1(x,y,0);

        % Calculate peak amplitude and full width at half max (FWHM)
        borderAmplitude = yPeak - yOuterEdge;
        yHalfMax = yOuterEdge + 0.5*borderAmplitude;
        xHalfMaxOuter = interp1(y(x < xPeak),x(x < xPeak),yHalfMax);
        xHalfMaxInner = x(find(x > xPeak & y <= yHalfMax,1,'first'));
        if isempty(xHalfMaxInner)
            xHalfMaxInner = NaN;
        end
        FWHM = xHalfMaxInner - xHalfMaxOuter;

        % Calculate patch radius and mean amplitude of center of the patch
        lawnMask = bwdist(lawn.labeled ~= j);
        lawnMask = -lawnMask + max(lawnMask,[],'all');
        lawn.edge = lawn.edge | (lawnMask <= max(x)*scale);
        lawnAreaPixels = sum(lawnMask <= max(x)*scale,'all'); % (pixels)
        lawnRadius = sqrt(lawnAreaPixels/pi)/scale; % (mm)
        meanAmplitude = mean(lawn.imageNormalized(lawnMask <= max(x)*scale),'all') - ...
            yOuterEdge;
        centerAmplitude = mean(lawn.imageNormalized(lawnMask <= max(x)*scale/2),'all') - ...
            yOuterEdge;
        borderCenterRatio = borderAmplitude/centerAmplitude;

        % Combine analyses in table
        circularity = lawn.circularity(j);
        analysis = table(xPeak,yPeak,xOuterEdge,yOuterEdge,...
            xHalfMaxOuter,xHalfMaxInner,yHalfMax,FWHM,lawnRadius,circularity,...
            borderAmplitude,meanAmplitude,centerAmplitude,borderCenterRatio);
        lawn.analysis = [lawn.analysis;analysis];
    end
end
end