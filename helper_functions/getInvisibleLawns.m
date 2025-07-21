function [lawn] = getInvisibleLawns(lawnFileName,firstFrame,arena,orientation,lawnDiameter,lawnSpacing,OD600)
% [lawns] = GETINVISIBLELAWNS(lawnFileName,firstFrame,arena,orientation,
%   lawnDiameter,lawnSpacing,OD600,camera)
%
%   GETINVISIBLELAWNS creates an estimate of the locations of bacterial 
%   lawns in the behavioral arena. GETINVISIBLELAWNS expects a lawnFileName
%   corresponding to a video where a dark object (e.g. black cardstock) is
%   passed between the lightsource and the experimental plate, creating a
%   "darkfield" image of the plate. The cardstock scatters the light
%   creating increased contrast of the bacterial lawns, enabling automatic
%   identfication of the lawn edges. GETINVISIBLELAWNS expects that the
%   video lawnFileName is not necessarily in the same position and
%   orientation as the experimental video itself and thus uses the arena 
%   and orientation structures to translate, rotate, and scale the lawn
%   estimates. If not lawnFileName is provided, GETINVISIBLELAWNS will use
%   the lawnDiameter and lawnSpacing to estimate the location of the
%   lawn(s) using an isometric grid pattern.
%
%   INPUTS:
%       - lawnFileName [str]:  file name of a video file (e.g. 
%           YYYY-MM-DD_HH-MM-SS_C.avi); C here is used to denote the camera
%           number; if no lawn contrast video exists, pass an empty array
%           (i.e. []) and the lawn locations will be estimated
%       - firstFrame [MxN uint8]: image of the experiment video;
%           often the first frame of the contrast video being analyzed
%       - arena [struct]: structure created by getAcetateArena() for
%           the video file containing information about the properties of
%           the acetate arena in the experiment video
%       - orientation [struct]: structure created by
%           getAcetateOrientation() for the video file containing
%           information about the orientation of the acetate arena in the
%           experiment video
%       - lawnDiameter [double]: the diameter of the lawn(s) in mm (e.g.
%           1.5)
%       - lawnSpacing [double]: the center-to-center distance between lawns
%           in mm (e.g. 6); use 0 for a single lawn
%       - OD600 [double]: the OD600 concentration of the lawn(s) (e.g. 0.5)
%       - camera [double]: camera identifier for which the experiment video
%           recording was taken (e.g. 1 or 2)
%
%   OUTPUTS:
%       - lawn [struct]: a structure containing the following fields
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
%   Written 5/2/2022 by Jess Haley in MATLAB R2022a.
%
%   See also ANALYZEFORAGING, GETVIDEOINFO, GETACETATEARENA,
%   GETACETATEORIENTATION, GETLAWNSTEMPLATE, GETLAWNSTHRESHOLD,
%   GETLAWNSFILTER, GETLAWNSCIRCLE, REGISTERIMAGE, REGIONPROPS, BWLABEL, 
%   IMDILATE, IMERODE, IMROTATE, IMTRANSLATE, PDIST.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Create Lawn Masks from First Frame

% If OD600 is a cell array (meaning multiple concentrations per plate)
if iscell(OD600)
    OD600_cell = OD600{:};
    OD600 = median(OD600{:});
end

% If OD600 is zero, make non-zero
OD600 = max(OD600,1e-10);

% Get and adjust lawn template
template = getLawnsTemplate(size(arena.mask),arena.scale,...
    arena.diameter,lawnDiameter,lawnSpacing);
lawn.templateB = transformLabeledImage(template,flip(arena.offset),...
    orientation.angle,1);
lawnRadius = mean(lawn.templateB.radii);
numLawn = length(lawn.templateB.radii);

if OD600 > 0.01

    % Use threshold method on behavior video
    lawn.thresholdB = getLawnsThreshold(firstFrame,arena.mask,lawnRadius);
    if ~isempty(lawn.thresholdB.centers)
        lawn.thresholdB.dist = min(pdist2(lawn.thresholdB.centers,lawn.templateB.centers),[],2);
    else
        lawn.thresholdB.dist = [];
    end

    % Update template to threshold fits or orientation dot
    if length(lawn.thresholdB.radii) > 2 & mean(lawn.thresholdB.dist) < arena.scale
        templateMask = registerImage(lawn.thresholdB.mask,lawn.templateB.mask,...
            lawn.thresholdB.centers,lawn.templateB.centers);
        lawn.templateB = getLawnProperties(templateMask);
    else
        lawn.templateB = transformLabeledImage(lawn.templateB,-arena.offset,orientation.angle,1);
    end

    % Update distances to new template
    if ~isempty(lawn.thresholdB.dist)
        lawn.thresholdB.dist = min(pdist2(lawn.thresholdB.centers,lawn.templateB.centers),[],2);
    end

    % Use circle method on behavior video
    lawn.circleB = getLawnsCircle(firstFrame,arena.mask,lawn.templateB.mask);
    if ~isempty(lawn.circleB.centers)
        lawn.circleB.dist = min(pdist2(lawn.circleB.centers,lawn.templateB.centers),[],2);
    else
        lawn.circleB.dist = [];
    end

    % Adjust template if not yet done
    if (length(lawn.thresholdB.radii) <= 2 | mean(lawn.thresholdB.dist) >= arena.scale) & numLawn > 1
        templateMask = registerImage(lawn.circleB.mask,lawn.templateB.mask,...
            lawn.circleB.centers,lawn.templateB.centers);
        lawn.templateB = getLawnProperties(templateMask);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Create Lawn Masks from Contrast Video

% If there is a lawn video file, try threshold, filter, and circle methods
if ~isempty(lawnFileName) && OD600 > 0.01
    
    % Get lawn video, arena, and orientation
    lawn.video = getVideoInfo(lawnFileName);
    lawn.arena = getAcetateArena(lawn.video.firstFrame,arena.diameter);
    if isnan(orientation.distance)
        lawn.orientation = orientation;
    else
        lawn.orientation = getAcetateOrientation(lawn.video.firstFrame,lawn.arena);
    end

    % Get and adjust lawn template
    lawn.template =  getLawnsTemplate(size(lawn.arena.mask),lawn.arena.scale,...
        lawn.arena.diameter,lawnDiameter,lawnSpacing);
    lawnRadius = mean(lawn.template.radii);
    numLawn = length(lawn.template.radii);

    % Get image registration transformation
    [~,lawn.registration] = registerImage(firstFrame,lawn.video.firstFrame);
    arenaReg = imwarp(lawn.arena.mask,lawn.registration,OutputView=imref2d(size(arena.mask)));
    lawn.regFit = sum(arenaReg & arena.mask,'all')/sum(arena.mask,'all');
    
    % Detect lawns using threshold from first frame of lawn video file
    lawn.threshold = getLawnsThreshold(lawn.video.firstFrame,lawn.arena.mask,...
        lawnRadius);

    % Detect lawns using filter from lawn video file
    lawn.filter = getLawnsFilter(lawn.video.vidObject,lawn.arena.mask,...
        lawnRadius,numLawn);
    
    % Get fit of threshold method by comparing distances to template lawns
    if ~isempty(lawn.threshold.centers)
        lawn.threshold.dist = min(pdist2(lawn.threshold.centers,lawn.template.centers),[],2);
    else
        lawn.threshold.dist = [];
    end

    % Get fit of filter method by comparing distances to template lawns
    if ~isempty(lawn.filter.centers)
        lawn.filter.dist = min(pdist2(lawn.filter.centers,lawn.template.centers),[],2);
    else
        lawn.filter.dist = [];
    end

    % Use best fit to transform template mask to better match found patches
    lawn.template = transformLabeledImage(lawn.template,flip(lawn.arena.offset),...
        lawn.orientation.angle,1);
    if (mean(lawn.threshold.dist) < mean(lawn.filter.dist) |...
            size(lawn.filter.centers,1) <= 2) & size(lawn.threshold.centers,1) > 2 && ...
            mean(lawn.threshold.dist) < lawn.arena.scale
        templateMask = registerImage(lawn.threshold.mask,lawn.template.mask,...
            lawn.threshold.centers,lawn.template.centers);
        lawn.template = getLawnProperties(templateMask);
    elseif size(lawn.filter.centers,1) > 2
        templateMask = registerImage(lawn.filter.mask,lawn.template.mask,...
            lawn.filter.centers,lawn.template.centers);
        lawn.template = getLawnProperties(templateMask);
    end

    % Find circles using updated template
    lawn.circle = getLawnsCircle(lawn.video.vidObject,lawn.arena.mask,lawn.template.mask);

    % Adjust template if not yet done
    if size(lawn.threshold.centers,1) <= 2 & size(lawn.filter.centers,1) <= 2 & ...
            size(lawn.circle.centers,1) > 2
        templateMask = registerImage(lawn.circle.mask,lawn.template.mask,...
            lawn.circle.centers,lawn.template.centers);
        lawn.template = getLawnProperties(templateMask);
    end

     % Update distances to new template
    if ~isempty(lawn.threshold.dist)
        lawn.threshold.dist = min(pdist2(lawn.threshold.centers,lawn.template.centers),[],2);
    end
    if ~isempty(lawn.filter.dist)
        lawn.filter.dist = min(pdist2(lawn.filter.centers,lawn.template.centers),[],2);
    end
    if ~isempty(lawn.circle.centers)
        lawn.circle.dist = min(pdist2(lawn.circle.centers,lawn.template.centers),[],2);
    else
        lawn.circle.dist = [];
    end
else
    lawn.regFit = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Try to use threshold and circles from behavior video

% Start with lawn template
if isfield(lawn,'template')
    lawnLabeled = imwarp(lawn.template.labeled,lawn.registration,'nearest',...
        Outputview=imref2d(size(firstFrame)));
else
    lawnLabeled = lawn.templateB.labeled;
end
lawnMethod = repmat({'template'},numLawn,1);

if OD600 > 0.1

    % Add thresholds if good fit
    goodFit = find(lawn.thresholdB.dist < arena.scale & lawn.thresholdB.circularity > 0.89);
    for numROI = 1:length(goodFit)
        centerROI = round(lawn.thresholdB.centers(goodFit(numROI),:));
        replaceID = lawnLabeled(centerROI(2),centerROI(1));
        circleID = lawn.circleB.labeled(centerROI(2),centerROI(1));
        if replaceID > 0
            lawnLabeled(lawnLabeled == replaceID) = 0;
            % lawnLabeled(lawn.thresholdB.labeled == goodFit(numROI)) = replaceID;
            % lawnMethod{replaceID} = 'thresholdB';
            overlap = sum(lawn.thresholdB.labeled == goodFit(numROI) & lawn.circleB.labeled == circleID,'all')/...
                sum(lawn.thresholdB.labeled == goodFit(numROI) | lawn.circleB.labeled == circleID,'all');
            if overlap > 0.75
                lawnLabeled(lawn.thresholdB.labeled == goodFit(numROI) | lawn.circleB.labeled == circleID) = replaceID;
                lawnMethod{replaceID} = 'thresholdBAdjusted';
            else
                lawnLabeled(lawn.thresholdB.labeled == goodFit(numROI)) = replaceID;
                lawnMethod{replaceID} = 'thresholdB';
            end
        end
    end

    % Add circles if good fit
    goodFit = find(lawn.circleB.dist < arena.scale & lawn.circleB.strength > 0.1 & ...
        lawn.circleB.radii >= 0.9*arena.scale*lawnDiameter/2);
    for numROI = 1:length(goodFit)
        % centerROI = round(lawn.circleB.centers(goodFit(numROI),:));
        % replaceID = lawnLabeled(centerROI(2),centerROI(1));
        lawnObjects = regionprops(lawnLabeled,'Centroid');
        [~,replaceID] = min(pdist2(lawn.circleB.centers(goodFit(numROI),:),vertcat(lawnObjects.Centroid)));
        if replaceID > 0 && strcmp(lawnMethod{replaceID},'template')
            lawnLabeled(lawnLabeled == replaceID) = 0;
            lawnLabeled(lawn.circleB.labeled == goodFit(numROI)) = replaceID;
            lawnMethod{replaceID} = 'circleB';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Try to use threshold, filter, and circles from contrast video

if isfield(lawn,'registration')

    % Add thresholds if good fit
    goodFit = find(lawn.threshold.dist < arena.scale & lawn.threshold.circularity > 0.95);
    regMask = imwarp(lawn.threshold.mask,lawn.registration,'nearest',OutputView=imref2d(size(firstFrame)));
    regInfo = getLawnProperties(regMask);
    for numROI = 1:length(goodFit)
        [~,ind] = min(pdist2(lawn.threshold.centers(goodFit(numROI),:),regInfo.centers));
        centerROI = round(regInfo.centers(ind,:));
        replaceID = lawnLabeled(centerROI(2),centerROI(1));
        if replaceID > 0 && strcmp(lawnMethod{replaceID},'template')
            lawnLabeled(lawnLabeled == replaceID) = 0;
            lawnLabeled(regInfo.labeled == ind) = replaceID;
            lawnMethod{replaceID} = 'threshold';
        end
    end

    % Add filters if good fit (override other fits)
    goodFit = find(lawn.filter.dist < arena.scale & lawn.filter.circularity > 0.8);
    regMask = imwarp(lawn.filter.mask,lawn.registration,'nearest',OutputView=imref2d(size(firstFrame)));
    regLabeled = bwlabel(imwarp(lawn.circle.mask,lawn.registration,'nearest',OutputView=imref2d(size(firstFrame))));
    regInfo = getLawnProperties(regMask);
    for numROI = 1:length(goodFit)
        [~,ind] = min(pdist2(lawn.filter.centers(goodFit(numROI),:),regInfo.centers));
        centerROI = round(regInfo.centers(ind,:));
        % replaceID = lawnLabeled(centerROI(2),centerROI(1));
        lawnObjects = regionprops(lawnLabeled,'Centroid','Circularity');
        filterObjects = regionprops(regInfo.labeled == ind,'Circularity');
        [~,replaceID] = min(pdist2(regInfo.centers(ind,:),vertcat(lawnObjects.Centroid)));
        circleID = regLabeled(centerROI(2),centerROI(1));
        if (filterObjects.Circularity > lawnObjects(replaceID).Circularity) || ...
                ~contains(lawnMethod{replaceID},'threshold')
            lawnLabeled(lawnLabeled == replaceID) = 0;
            if circleID > 0 & numLawn > 1
                overlap = sum(regInfo.labeled == ind & regLabeled == circleID,'all')/...
                    sum(regInfo.labeled == ind | regLabeled == circleID,'all');
                if overlap > 0.75 || lawn.circle.strength(circleID) > 0.15
                    lawnLabeled(regInfo.labeled == ind | regLabeled == circleID) = replaceID;
                    lawnMethod{replaceID} = 'filterAdjusted';
                else
                    lawnLabeled(regInfo.labeled == ind) = replaceID;
                    lawnMethod{replaceID} = 'filter';
                end
            else
                lawnLabeled(regInfo.labeled == ind) = replaceID;
                lawnMethod{replaceID} = 'filter';
            end
        end
    end

    % Add circles if good fit
    goodFit = find(lawn.circle.dist < 0.8*arena.scale & lawn.circle.strength > 0.12);
    regMask = imwarp(lawn.circle.mask,lawn.registration,'nearest',OutputView=imref2d(size(firstFrame)));
    regInfo = getLawnProperties(regMask);
    for numROI = 1:length(goodFit)
        [~,ind] = min(pdist2(lawn.circle.centers(goodFit(numROI),:),regInfo.centers));
        lawnObjects = regionprops(lawnLabeled,'Centroid');
        [~,replaceID] = min(pdist2(regInfo.centers(ind,:),vertcat(lawnObjects.Centroid)));
        % centerROI = round(regInfo.centers(ind,:));
        % replaceID = lawnLabeled(centerROI(2),centerROI(1));
        if replaceID > 0 && strcmp(lawnMethod{replaceID},'template')
            lawnLabeled(lawnLabeled == replaceID) = 0;
            lawnLabeled(regInfo.labeled == ind) = replaceID;
            lawnMethod{replaceID} = 'circle';
        end
    end

    % Add circles with adjusted radius if okay fit
    okayFit = find(lawn.circle.dist < 0.8*arena.scale & lawn.circle.strength <= 0.12 & ...
        lawn.circle.strength > 0.05);
    for numROI = 1:length(okayFit)
        [~,ind] = min(pdist2(lawn.circle.centers(okayFit(numROI),:),regInfo.centers));
        lawnObjects = regionprops(lawnLabeled,'Centroid');
        [~,replaceID] = min(pdist2(regInfo.centers(ind,:),vertcat(lawnObjects.Centroid)));
        % centerROI = round(regInfo.centers(ind,:));
        % replaceID = lawnLabeled(centerROI(2),centerROI(1));
        if replaceID > 0 && strcmp(lawnMethod{replaceID},'template')
            lawnLabeled(lawnLabeled == replaceID) = 0;
            radiusDiff = mean(lawn.circle.radii(goodFit)) - lawn.circle.radii(okayFit(numROI));
            if radiusDiff < 1 || isempty(goodFit)
                lawnLabeled(regInfo.labeled == ind) = replaceID;
                lawnMethod{replaceID} = 'circle';
            else
                lawnMask = bwdist(regInfo.labeled == ind) <= radiusDiff;
                lawnLabeled(lawnMask) = replaceID;
                lawnMethod{replaceID} = 'circleAdjusted';
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. For remaining lawns, adjust template radii

lawnObjects = regionprops(lawnLabeled);
lawnRadii = sqrt(vertcat(lawnObjects.Area)./pi);
radiusUpdate =  mean(lawnRadii(~strcmp(lawnMethod,'template')));
for numROI = 1:numLawn
    if strcmp(lawnMethod{numROI},'template')
        radiusDiff = radiusUpdate - lawnRadii(numROI);
        if radiusDiff >= 1
            lawnMask = bwdist(lawnLabeled == numROI) <= radiusDiff;
            lawnLabeled(lawnMask) = numROI;
            lawnMethod{numROI} = 'templateAdjusted';
        end
    end
end

% Edit manual adjustment variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Query user to verify mask or manually adjust

% if numLawn == 1 && OD600 > 0.1
%     figure(2)
%     set(gcf,'Position',[700 700 350 350])
%     imshow(lawn.mask)
%     xlim(size(lawn.mask,2)/2 + [-250 250])
%     ylim(size(lawn.mask,1)/2 + [-250 250])
% 
%     figure(1)
%     set(gcf,'Position',[1000 100 1000 1000]);
%     imagesc(filter.image)
%     title('Press any key to confirm or skip changes')
%     xlim(size(lawn.mask,2)/2 + [-250 250])
%     ylim(size(lawn.mask,1)/2 + [-250 250])
%     translatedMask = imtranslate(lawn.mask,-arena.offset,'FillValues',0) > 0;
%     lawnObjects = regionprops(translatedMask,'Centroid','Area');
%     roi = drawcircle('Center',lawnObjects.Centroid,...
%         'Radius',sqrt(lawnObjects.Area/pi));
%     changeMask = 0;
%     w = NaN;
%     while w ~= 1
%         w = waitforbuttonpress;
%         if w == 0
%             changeMask = 1;
%         end
%     end
%     if changeMask == 1
%         manualMask = zeros(size(lawn.mask));
%         lawnCenters = round(roi.Center);
%         manualMask(lawnCenters(2),lawnCenters(1)) = 1;
%         manualMask = bwdist(manualMask) <= roi.Radius;
%         manualMask = imtranslate(manualMask,arena.offset,'FillValues',0) > 0;
%         lawn.mask = manualMask;
%         lawn.method = 'manual';
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Recenter lawn on worm video

lawn.mask = lawnLabeled > 0;
lawn.labeled = bwlabel(lawn.mask);

% Create mask labeling closest lawn to each pixel
closestROI = nan(size(lawn.labeled,1),size(lawn.labeled,2),max(lawn.labeled(:)));
for numROI = 1:max(lawn.labeled(:))
    closestROI(:,:,numROI) = bwdist(lawn.labeled == numROI);
end
[~,lawn.closest] = min(closestROI,[],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Get Lawn Info

% Get lawn image and mask
if exist('OD600_cell','var')
    lawn.OD600 = lawn.labeled;
    lawn.closestOD600 = lawn.closest;
    % lawn.nearLawnOD600 = lawn.nearLawnLabeled;
    for numROI = 1:max(lawn.labeled(:))
        templateCenter = round(template.centers(numROI,:));
        lawnID = lawn.closest(templateCenter(2),templateCenter(1));
        lawn.OD600(lawn.labeled == lawnID) = OD600_cell(numROI);
        lawn.closestOD600(lawn.closest == lawnID) = OD600_cell(numROI);
    end
    lawn.mask = lawn.OD600 > 0;
    lawn.labeled = bwlabel(lawn.mask);
    
    % Recreate mask labeling closest lawn to each pixel
    closestROI = nan(size(lawn.labeled,1),size(lawn.labeled,2),max(lawn.labeled(:)));
    for numROI = 1:max(lawn.labeled(:))
        closestROI(:,:,numROI) = bwdist(lawn.labeled == numROI);
    end
    [~,lawn.closest] = min(closestROI,[],3);
else
    lawn.OD600 = double(lawn.mask).*OD600;
    lawn.closestOD600 = ones(size(lawn.mask))*OD600;
end

% Get updated lawn centers and radii
lawnObjects = regionprops(lawn.labeled,'Centroid','Area','Circularity');
lawn.centers = vertcat(lawnObjects.Centroid);
lawn.radii = sqrt(vertcat(lawnObjects.Area)/pi);
lawn.circularity = vertcat(lawnObjects.Circularity);

% Add lawn method
lawn.method = cell(size(lawn.radii));
for numROI = 1:size(lawn.centers,1)
    labeledCenter = round(lawn.centers(numROI,:));
    lawn.method{numROI} = lawnMethod{lawnLabeled(labeledCenter(2),labeledCenter(1))};
end

% Get lawn statistics
[X,Y] = find(lawn.mask);
lawn.area = length(X)/(arena.scale^2); % mm^2
lawn.diameterMean = mean(2.*lawn.radii)/arena.scale; % mm
lawn.diameterSTD = std(2.*lawn.radii)/arena.scale; % mm
distances = pdist(lawn.centers)/arena.scale;
distances = distances(distances < 1.5*lawnSpacing);
lawn.spacingMean = mean(distances); % mm
lawn.spacingSTD = std(distances); % mm
lawn.distributionCenter = [mean(X),mean(Y)]; % pixels
lawn.distributionSTD = std(sqrt((X - mean(X)).^2 + ...
    (Y - mean(Y)).^2))/arena.scale; % mm

end