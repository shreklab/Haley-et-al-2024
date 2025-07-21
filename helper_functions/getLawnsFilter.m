function [filter] = getLawnsFilter(vidObject,arenaMask,lawnRadius,numLawn)
% [filter] = GETLAWNSFILTER(image,arenaMask,lawnRadius,numLawn)
%
%   GETLAWNSFILTER creates an estimate of the locations of bacterial 
%   lawns in the behavioral arena given a VideoReader object (ideally 
%   for a video without worms) and information about the expected position 
%   of the lawns using image processing and morphological operations. 
%   GETLAWNSFILTER expects a video where a dark object (e.g. black 
%   cardstock) is passed between the lightsource and the experimental plate,
%   creating a "darkfield" image of plate. The cardstock scatters the light
%   creating increased contrast of the bacteral lawns, enabling automatic
%   identification of the lawn edges. GETLAWNSFILTER expects that the video
%   used contains an arena that is not necessarily in the same position and
%   orientation as the experimental video itself and enables translation, 
%   rotation, and scaling of the lawn mask to match the experimental video.
%
%   INPUTS:
%       - vidObject [1x1 VideoReader]: video object used with read() to
%           get frames of the contrast video
%       - arenaMask [1024x1024 logical]: mask of the arena (1 = inside 
%           arena); if no arena, use false(size(image))
%       - lawnRadius [double]: radii of the expected lawn(s) in pixels 
%           (e.g. 24.5)
%       - numLawn [double]: number of lawns expected in the arena (e.g. 19)
%       - translate [1x2 double]: the number of pixels in X,Y that the 
%           image needs to be translated by in order to center the arena
%           (e.g. [20.7, -4.2])
%       - rotate [double]: the number of degrees the arena arena needs to
%           be rotated by in order to center the reference orientation mark
%           in its expected position (e.g. -3.1)
%       - resize [double]: the scale at which the image needs to be resized
%           by (e.g. 1.2); applicable only when 2 cameras have been used in 
%           the same experiment
%
%   OUTPUTS:
%       - filter [struct]: a structure containing the following fields
%           - mask [1024x1024  logical]: mask of the lawn
%               locations (1 = lawn)
%           - labeled [1024x1024 double]: image of the lawn locations with 
%               each lawn uniquely labeled with a number in the set 
%               [1,2,...,numLawn]
%           - image [1024x1024 double]: maximum intensity projection of the 
%               filtered image across all frames of vidObject
%           - centers [numLawnx2 double]: X,Y positions of each lawn in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - radii [numLawnx1 double]: radii of each lawn in pixels with
%               row index corresponding to the label(s) in labeled
%
%   Written 4/13/2022 by Jess Haley in MATLAB R2021a.
%
%   See also ANALYZEFORAGING, GETINVISIBLELAWNS, VIDEOREADER, BWLABEL, 
%   REGIONPROPS, IMBINARIZE, IMFILTER, IMERODE, IMFILL, IMTRANSLATE, 
%   IMRESIZE, IMCROP, IMROTATE, PADARRAY, CONV2, STREL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(vidObject)

    % Start parallel pool
    % if isempty(gcp('nocreate'))
    %     pool = parpool(24);
    % end
   
    % Get expected lawn area
    lawnArea = pi*lawnRadius^2;
    
    % Get mask of visible lawns for each frame
    lawnMask = zeros([vidObject.Width,vidObject.Height,vidObject.NumFrames]);
    lawnImage = zeros([vidObject.Width,vidObject.Height,vidObject.NumFrames]);
    for k = 1:vidObject.NumFrames
        % Get current frame
        currentFrame = double(read(vidObject,k));
        
        % Fitler and binarize frame
        filteredFrame = imfilter(currentFrame,fspecial('average',4)) - ...
            imfilter(currentFrame,fspecial('average',40));
        binaryFrame = imbinarize(filteredFrame);

        % Get lawn objects
        if numLawn == 1
            labeledFrame = bwlabel(binaryFrame.*~(bwdist(~arenaMask) < 30));
            frameObjects = regionprops(labeledFrame,'FilledArea','Eccentricity');
            ind = find([frameObjects.FilledArea] > 1000 & ...
                [frameObjects.Eccentricity] < 0.85);
            lawnReal = false(size(currentFrame));
            for i = 1:length(ind)
                lawnReal(labeledFrame == ind(i)) = true;
            end
        else
            % Get darkfield image
            darkfieldFrame = binaryFrame.*~(bwdist(~(currentFrame < 100)) < 10).*...
                ~(bwdist(~arenaMask) < 30); % where the paper was
            labeledFrame = bwlabel(darkfieldFrame);
            frameObjects = regionprops(labeledFrame,'ConvexArea','Eccentricity');
            ind = find([frameObjects.ConvexArea] > lawnArea*0.5 & ...
                [frameObjects.Eccentricity] < 0.85);
            lawnReal = false(size(currentFrame));
            for i = 1:length(ind)
                lawnReal(labeledFrame == ind(i)) = true;
            end
            
            % Get brightfield image
            brightfieldFrame = ~(bwdist(binaryFrame)<=1).*...
                ~(bwdist(~(currentFrame > 100)) < 10).*...
                ~(bwdist(~arenaMask) < 30);
            labeledFrame = bwlabel(brightfieldFrame);
            frameObjects = regionprops(labeledFrame,'MajorAxis','Area',...
                'Eccentricity','ConvexArea');
            ind = find([frameObjects.Area] > 1000 & ...
                [frameObjects.MajorAxisLength] < 2*lawnRadius*2 & ...
                [frameObjects.Eccentricity] < 0.85 & ...
                [frameObjects.ConvexArea] > lawnArea*0.7);
            for i = 1:length(ind)
                lawnReal(labeledFrame == ind(i)) = true;
            end
        end
        lawnMask(:,:,k) = lawnReal;
        lawnImage(:,:,k) = filteredFrame;
    end
    
    % Sum masked image objects across frames and fill holes
    summedImage = imfill(sum(lawnMask,3) > 0,'holes');

    % Check if mask is circular or needs to be connected
    % if numLawn == 1
    %     if sum(summedImage(:)) < 0.7*lawnArea
    %         [center,radius] = imfindcircles(summedImage,round(lawnRadius*[0.7,1.3]),...
    %             'Sensitivity',0.95);
    %         if length(radius) == 1
    %             X = repmat(1:size(arenaMask,1),size(arenaMask,2),1)';
    %             Y = repmat(1:size(arenaMask,2),size(arenaMask,1),1);
    %             distances = sqrt((X-center(2)).^2 + ...
    %                 (Y-center(1)).^2);
    %             summedImage(distances <= radius) = 1;
    %         end
    %     end
    % end
    
    % Check if mask has lawn object(s) that are circular and correct size
    labeledImage = bwlabel(summedImage);
    imageObjects = regionprops(labeledImage,'Area','Circularity','Eccentricity');
    objectRadii = sqrt(vertcat(imageObjects.Area)/pi);
    objectCircularity = vertcat(imageObjects.Circularity);
    objectEccentricity = vertcat(imageObjects.Eccentricity);
    lawnObjects = find(objectRadii > 0.6*lawnRadius & ...
        objectRadii < 2*lawnRadius & ...
        objectCircularity > 0.2 & objectEccentricity < 0.85);
    summedImage = false(size(summedImage));
    for i = 1:length(lawnObjects)
        summedImage(labeledImage == lawnObjects(i)) = 1;
    end

    % Remove nubs
    nublessImage = bwdist(bwdist(~summedImage)>=10)<=10;
    summedImage = (summedImage - nublessImage == 0) & (summedImage > 0);
    
    % Blur image via convolution to smooth masks
    windowSize = 10;
    kernel = ones(windowSize)/windowSize^2;
    blurryImage = conv2(single(summedImage),kernel,'same');
    smoothedImage = blurryImage > 0.5;

    % Check objects again
    labeledImage = bwlabel(smoothedImage);
    imageObjects = regionprops(labeledImage,'Area');
    lawnObjects = find(sqrt([imageObjects.Area]/pi) > 0.6*lawnRadius);
    smoothedImage = false(size(smoothedImage));
    for i = 1:length(lawnObjects)
        smoothedImage(labeledImage == lawnObjects(i)) = 1;
    end
    
    mask = smoothedImage;
else
    mask = false(size(arenaMask));
    lawnImage = zeros([size(arenaMask),2]);
end

% Label lawn objects and get centers and radii
filter = getLawnProperties(mask,max(lawnImage,[],3));

end

