function [circle] = getLawnsCircle(vidObject,arenaMask,templateMask)
% [circle] = GETLAWNSCIRCLE(vidObject,arenaMask,templateMask)
%
%   GETLAWNSCIRCLE creates an estimate of the locations of bacterial 
%   lawns in the behavioral arena given a VideoReader object (ideally 
%   for a video without worms) and information about the expected position 
%   of the lawns using image processing and morphological operations.
%   GETLAWNSCIRCLE expects a video where a dark object (e.g. black 
%   cardstock) is passed between the lightsource and the experimental plate,
%   creating a "darkfield" image of plate. The cardstock scatters the light
%   creating increased contrast of the bacteral lawns, enabling automatic
%   identification of the lawn edges. GETLAWNSCIRCLE uses IMFINDCIRCLES to 
%   identify the putative lawns on every frame and creates a weighted 
%   average image. GETLAWNSCIRCLE expects that the video used contains an 
%   arena that is not necessarily in the same position and orientation as 
%   the experimental video itself and enables translation, rotation, and 
%   scaling of the lawn mask to match the experimental video.
%
%   INPUTS:
%       - vidObject [1x1 VideoReader]: video object used with read() to
%           get frames of the contrast video
%       - arenaMask [1024x1024 logical]: mask of the arena (1 = inside 
%           arena); if no arena, use false(size(image))
%       - templateMask [1024x1024 logical]: mask of the expected lawn
%           locations (1 = lawn)
%
%   OUTPUTS:
%       - circle [struct]: a structure containing the following fields
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
%           - strength [numLawnx1 double]: the maximum relative strength of
%               each lawn as defined by IMFINDCIRCLES
%
%   Written 2/14/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEFORAGING, GETINVISIBLELAWNS, GETLAWNSTEMPLATE,
%   IMFINDCIRCLES, VIDEOREADER, BWLABEL, REGIONPROPS, IMBINARIZE, IMFILTER,
%   IMERODE, IMFILL, IMTRANSLATE, IMRESIZE, IMCROP, IMROTATE, PADARRAY, 
%   IMBINARIZE, CONV2, STREL, SUB2IND, BWDIST, SQUAREFORM, MESHGRID.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(vidObject)

    if isa(vidObject,'VideoReader')
        lastFrame = nan(vidObject.Width,vidObject.Height);
    else
        lastFrame = double(vidObject);
        vidObject = struct('NumFrames',1,'Width',size(lastFrame,1),...
            'Height',size(lastFrame,2));
    end

    templateLabeled = bwlabel(templateMask);
    templateDilate = bwdist(templateMask) <= 0.25*max(bwdist(~templateMask),[],'all');
    numLawn = sum(unique(templateLabeled) > 0);
    lawnObjects = regionprops(templateLabeled,'Area');
    lawnRadius = sqrt(mean([lawnObjects(:).Area])/pi);

    % Start parallel pool
    % if isempty(gcp('nocreate'))
    %     pool = parpool(24);
    % end

    % Get mask of visible lawns for each frame
    lawnCircles = cell(vidObject.NumFrames,1);
    lawnImage = zeros([vidObject.Width,vidObject.Height,vidObject.NumFrames]);
    % currentFrame = nan(vidObject.Width,vidObject.Height);
    % figure
    for k = 1:vidObject.NumFrames

        % Get current frame
        if isa(vidObject,'VideoReader')
            currentFrame = double(read(vidObject,k));
        else
            currentFrame = lastFrame;
        end

        % Fitler and binarize frame
        filteredFrame = imfilter(currentFrame,fspecial('average',4)) - ...
            imfilter(currentFrame,fspecial('average',40));

        % Find circles in brightfield image
        brightfieldFrame = filteredFrame.*~(bwdist(~arenaMask) < 30).*...
            ~(bwdist(~(currentFrame > 100)) < 30);
        [centerB,radiusB,strengthB] = imfindcircles(brightfieldFrame,...
            round(lawnRadius*[0.7,1.7]),'Sensitivity',0.95,'ObjectPolarity','Dark');

        % Find circles in darkfield image
        darkfieldFrame = currentFrame.* ~(bwdist(~arenaMask) < 30).*...
            ~(bwdist(~(currentFrame < 100)) < 10);
        [centerD,radiusD,strengthD] = imfindcircles(darkfieldFrame,...
            round(lawnRadius*[0.7,1.7]),'Sensitivity',0.95,'ObjectPolarity','Bright');

        % Merge circles
        center = [centerB;centerD];
        radius = [radiusB;radiusD];
        strength = [strengthB;strengthD];

        % If still frame used and little to no lawns found, try to binarize
        if ~isa(vidObject,'VideoReader') & length(radius)/numLawn < 0.4
            binaryFrame = imbinarize(brightfieldFrame,'adaptive','Sensitivity',0.6);
            [center,radius,strength] = imfindcircles(binaryFrame,...
                round(lawnRadius*[0.7,1.7]),'Sensitivity',0.95,'ObjectPolarity','Dark');
        end

        % If only one lawn and none found, try to use filtered frame
        if numLawn == 1 & isempty(radius)
            binaryFrame = imbinarize(filteredFrame);
            [center,radius,strength] = imfindcircles(binaryFrame,...
                round(lawnRadius*[0.7,1.7]),'Sensitivity',0.95,'ObjectPolarity','Bright');
        end

        % Remove circles that don't overlap with template
        if ~isempty(center)
            [ind] = sub2ind(size(currentFrame),round(center(:,2)),round(center(:,1)));
            indOverlap = templateDilate(ind);
            center(~indOverlap,:) = []; radius(~indOverlap) = []; strength(~indOverlap) = [];
        end

        lawnCircles{k,:} = [center,radius,strength];
        lawnImage(:,:,k) = filteredFrame;
    end

    plotCircles = vertcat(lawnCircles{:});
    % figure; imshowpair(max(lawnImage,[],3),templateDilate); hold on;
    % for i = 1:size(plotCircles,1)
    %     viscircles(plotCircles(i,1:2),plotCircles(i,3),'LineWidth',3*plotCircles(i,4)/max(plotCircles(:,4)));
    % end

    if size(plotCircles,1) > 0
        distCircles = squareform(pdist(plotCircles(:,1:2)));
        indOverlap = [];
        if ~isempty(distCircles)
            [indOverlap(:,1),indOverlap(:,2)] = ...
                find(distCircles < 4*plotCircles(:,3) & distCircles > 0);
            indOverlap = sort(indOverlap,2);
        end
    noOverlap = ~ismember(1:size(plotCircles,1),unique(indOverlap(:)));
    bestCenter = plotCircles(noOverlap,1:2);
    bestRadius = plotCircles(noOverlap,3);
    bestStrength = plotCircles(noOverlap,4);
    while ~isempty(indOverlap)
        checkOverlap = indOverlap(indOverlap(:,1) == indOverlap(1,1),:);
        checkOverlap = unique(checkOverlap(:));

        % Weighted mean of center and radius by strength of circle fit
        % [maxStrength,bestFit] = max(plotCircles(checkOverlap,4));
        % bestCenter = [bestCenter;plotCircles(checkOverlap(bestFit),1:2)];
        % bestRadius = [bestRadius;plotCircles(checkOverlap(bestFit),3)];
        % rangeStrength = prctile(plotCircles(checkOverlap,4),[0 100]);
        % bestFit = find(plotCircles(checkOverlap,4) > rangeStrength(2)-0.5*diff(rangeStrength));
        maxStrength = prctile(plotCircles(checkOverlap,4),100);
        maxRadius = prctile(plotCircles(checkOverlap,3),100);
        bestFit = find(plotCircles(checkOverlap,4) > 0.8*maxStrength & ...
            plotCircles(checkOverlap,3) > 0.1*maxRadius);
        % bestFit = find(plotCircles(checkOverlap,4) > 0.05);
        [~,ind] = max(plotCircles(checkOverlap(bestFit),3));
        bestFit = bestFit(ind);
        bestCenter = [bestCenter;...
            sum(plotCircles(checkOverlap(bestFit),1:2).*plotCircles(checkOverlap(bestFit),4),1)./...
            sum(plotCircles(checkOverlap(bestFit),4),1)];
        bestRadius = [bestRadius;...
            sum(plotCircles(checkOverlap(bestFit),3).*plotCircles(checkOverlap(bestFit),4))./...
            sum(plotCircles(checkOverlap(bestFit),4))];
        bestStrength = [bestStrength;max(plotCircles(checkOverlap(bestFit),4))];
        
        plotCircles(checkOverlap,:) = [];
        if isempty(plotCircles) || size(plotCircles,1) == 1
            indOverlap = [];
            continue
        end

        distCircles = squareform(pdist(plotCircles(:,1:2)));
        indOverlap = [];
        [indOverlap(:,1),indOverlap(:,2)] = ...
            find(distCircles < 4*plotCircles(:,3) & distCircles > 0);
        indOverlap = sort(indOverlap,2);
    end

    distCircles = squareform(pdist(bestCenter));
    indOverlap = [];
    [indOverlap(:,1),indOverlap(:,2)] = find(distCircles < 2*bestRadius & distCircles > 0);
    if all(size(indOverlap) == [2,2])
        indOverlap = unique(indOverlap);
        [~,ind] = min(bestStrength(indOverlap));
        bestCenter(indOverlap(ind),:) = [];
        bestRadius(indOverlap(ind)) = [];
        bestStrength(indOverlap(ind)) = [];
    end

    % figure; imagesc(max(lawnImage,[],3));hold on;
    % viscircles(bestCenter,bestRadius);

    lawnMask = false(size(arenaMask));
    [X,Y] = meshgrid(1:size(arenaMask,1),1:size(arenaMask,2));
    for numROI = 1:length(bestRadius)
        distances = sqrt((X-bestCenter(numROI,1)).^2 + ...
            (Y-bestCenter(numROI,2)).^2);
        lawnMask(distances <= bestRadius(numROI)) = 1;
    end
    % figure; imshowpair(lawnMask,templateMask)

    lawnLabeled = bwlabel(lawnMask);
    lawnKeep = nan(numLawn,1);
    for numROI = 1:numLawn
        overlap = lawnLabeled.*(templateLabeled == numROI);
        [~,lawnKeep(numROI)] = max(histcounts(overlap,[1:length(bestRadius)+1]-0.5));
    end
    smoothedImage = false(size(lawnMask));
    for numROI = 1:numLawn
        smoothedImage(lawnLabeled == lawnKeep(numROI)) = 1;
    end
    % figure; imshow(smoothedImage)
    else
        smoothedImage = false(size(arenaMask));
        bestCenter = [];
        bestStrength = [];
    end
    circle.mask = smoothedImage;
else
    circle.mask = false(size(arenaMask));
    lawnImage = zeros([size(arenaMask),2]);
end

% Label lawn objects and get centers and radii
circle.labeled = bwlabel(circle.mask);
circle.image = max(lawnImage,[],3);
imageObjects = regionprops(circle.labeled,'Centroid','Area');
circle.centers = vertcat(imageObjects.Centroid);
circle.radii = sqrt(vertcat(imageObjects.Area)/pi);

% Get strength of fit
[~,ind] = min(pdist2(bestCenter,circle.centers),[],2);
circle.strength = bestStrength(ind);

end

