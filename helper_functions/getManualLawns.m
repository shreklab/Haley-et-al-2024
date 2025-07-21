function [lawn] = getManualLawns(lawnFile,lawn,scale)
% [lawn] = GETMANUALLAWNS(lawnFile,lawn,scale)
%
%   GETLAWNSMANUAL loads an image file containing the mask of the lawn 
%   locations and updates the lawn properties to correspond with the image.
%
%   INPUTS:
%       - lawnFile [str]: file name of an image file (e.g., .png)
%       - lawn [struct]: a structure output by GETMANUALLAWNS or
%           GETINVISIBLELAWNS (see below).
%       - scale [double]: calculated scale of the image in pixels/mm
%
%   OUTPUTS:
%        lawn [struct]: a structure containing the following fields
%           - video [struct]: structure created using getVideoInfo() on the
%               lawn file; contains metadata information about the video
%           - arena [struct]: structure created using getAcetateArena() on
%               the lawn file; contains metadata information about the
%               arena
%           - orientation [struct]: structure created using
%               getAcetateOrientation(); contains metadata about the
%               orientation of the arena
%           - template [struct]: structure created using
%               getLawnsTemplate(); contains a mask and its properties for
%               an estimate of lawn location(s) given expected spacing and
%               diameter
%           - threshold [struct]: structure created using
%               getLawnsThreshold(); contains a mask and its properties for
%               the lawn loacation(s) by using imaging processing
%               techniques on the first frame of the lawn file
%           - filter [struct]: structure created using
%               getLawnsFilter(); contains a mask and its properties for
%               the lawn loacation(s) by using imaging processing
%               techniques on all frames of the contrast video in lawn file
%           - mask [MxN logical]: mask of the lawn locations (1 = lawn)
%           - labeled [MxN double]: image of the lawn locations with 
%               each lawn uniquely labeled with a number in the set
%               [1,2,...,numLawn]
%           - method [str]: identifier for method used to estimate lawns
%               (i.e. 'template', 'threshold', or 'filter'); method depends
%               on existance and quality of lawn contrast video
%           - centers [numLawnx2 double]: X,Y positions of each lawn in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - radii [numLawnx1 double]: radii of each lawn in pixels with
%               row index corresponding to the label(s) in labeled
%           - closest [MxN double]: image showing the identity of the 
%               nearest lawn [1,2,...,numLawn] to each pixel
%           - circularity [numLawnx1 double]: circularity of each lawn with
%               row index corresponding to the label(s) in labeled (1 =
%               perfect circle)
%           - OD600 [MxN double]: image of the lawn(s) with pixel
%               value equal to the OD600 concentration (0 = off lawn)
%           - image [MxN uint8]: image of the lawn(s) with pixel
%               value scaled by concentration; image displays nicely for 
%               OD600 values between 0.05 and 10 using imshow()
%           - area [double]: summed area of lawn(s) in mm^2 (e.g. 40.2)
%           - diamaterMean [double]: mean diameter of the lawn(s) in mm
%               (e.g. 1.65)
%           - diameterSTD [double]: standard deviation of the diameter of
%               the lawn(s) in mm (e.g. 0.1)
%           - spacingMean [double]: mean center-to-center spacing of all
%               lawn(s) in mm (e.g. 6.01)
%           - spacingSTD [double]: standard deviation of the
%               center-to-center spacing of all lawn(s) in mm (e.g. 0.05)
%           - distributionCenter [1x2 double]: mean X,Y of all lawn(s) in
%               pixels (e.g. [512.0,512.0])
%           - distributionSTD [double]: standard devation of pixel
%               distances from all lawn(s) to the distributionCenter in mm 
%               (e.g. 3.2)
%
%   Written 3/4/2024 by Jess Haley in MATLAB R2024a.
%
%   See also ANALYZEFORAGING, GETINVISIBLELAWNS, IMREAD, BWLABEL,
%   REGIONPROPS, BWDIST.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy old masks
oldMask = lawn.mask;
oldLabeled = lawn.labeled;
oldMethod = lawn.method;

% Update masks
manualMask = imread(lawnFile);
lawn.mask = manualMask ~= 1; % Or ==0?
lawn.labeled = bwlabel(lawn.mask);

% Remove any extranneous pixels
lawnObjects = regionprops(lawn.labeled,'Area');
extraPixels = find(vertcat(lawnObjects.Area) < 10);
if ~isempty(extraPixels)
    for i = 1:length(extraPixels)
        lawn.mask(lawn.labeled == extraPixels(i)) = false;
    end
    lawn.labeled = bwlabel(lawn.mask);
end

% Create mask labeling closest lawn to each pixel
closestROI = nan(size(lawn.labeled,1),size(lawn.labeled,2),max(lawn.labeled(:)));
for numROI = 1:max(lawn.labeled(:))
    closestROI(:,:,numROI) = bwdist(lawn.labeled == numROI);
end
[~,lawn.closest] = min(closestROI,[],3);

% Get lawn OD600
lawn.OD600 = lawn.closestOD600.*lawn.mask;
for numROI = 1:max(lawn.labeled(:))
    lawn.closestOD600(lawn.closest == numROI) = unique(lawn.OD600(lawn.labeled == numROI));
end

% Get updated lawn centers and radii
lawnObjects = regionprops(lawn.labeled,'Centroid','Area','Circularity');
lawn.centers = vertcat(lawnObjects.Centroid);
lawn.radii = sqrt(vertcat(lawnObjects.Area)/pi);
lawn.circularity = vertcat(lawnObjects.Circularity);

% Add lawn method
for numROI = 1:max(lawn.labeled(:))
    centerROI = round(lawn.centers(numROI,:));
    oldID = oldLabeled(centerROI(2),centerROI(1));
    overlap = sum(oldLabeled == oldID & lawn.labeled == numROI,'all')/...
        sum(oldLabeled == oldID | lawn.labeled == numROI,'all');
    if overlap < 0.999
        lawn.method{numROI} = 'manual';
    else
        lawn.method{numROI} = oldMethod{oldID};
    end
end

% Get lawn statistics
[X,Y] = find(lawn.mask);
lawn.area = length(X)/(scale^2); % mm^2
lawn.diameterMean = mean(2.*lawn.radii)/scale; % mm
lawn.diameterSTD = std(2.*lawn.radii)/scale; % mm
distances = pdist(lawn.centers)/scale;
% distances = distances(distances < 1.5*lawnSpacing);
lawn.spacingMean = mean(distances); % mm
lawn.spacingSTD = std(distances); % mm
lawn.distributionCenter = [mean(X),mean(Y)]; % pixels
lawn.distributionSTD = std(sqrt((X - mean(X)).^2 + ...
    (Y - mean(Y)).^2))/scale; % mm

end