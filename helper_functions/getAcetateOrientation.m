function [orientation] = getAcetateOrientation(image,arena)
% [orientation] = GETACETATEORIENTATION(image,arena)
%
%   GETACETATEORIENTATION takes in an image containing a piece of acetate 
%   with both a large hole (i.e. an arena) to contain worm behavior and a
%   small hole cut to orient the pipetting template so that the location of
%   "invisible" (i.e. dilute E.coli) lawns can be estimated. 
%   GETACETATEORIENTATION returns a mask of the orientation mark and its 
%   properties using image processing and morphological operations.
%
%   INPUTS:
%       - image [1024x1024 uint8]: image of the acetate arena; often the
%           first frame of the video being analyzed
%       - arena [struct]: structure containing the mask and properties of
%           the acetate arena
%
%   OUTPUTS:
%       - orientation [struct]: a structure containing the following fields
%           - mask [1024x1024 logical]: mask of the reference dot
%               (1 = inside dot)
%           - center [1x2 double]: X,Y location of the reference
%               dot centroid in pixels
%           - distance [double]: Euclidean distance between the 
%               center of the arena and the center of the reference dot in 
%               pixels
%           - angle [double]: angle offset between the expected reference
%               dot position and the measured position given in degrees
%    
%   Written 4/13/2022 by Jess Haley in MATLAB R2021a.
%
%   See also ANALYZEFORAGING, GETACETATEARENA, REGIONPROPS, BWLABEL,
%   IMBINARIZE, IMDILATE, IMFILL, IMCOMPLEMENT, IMFILTER, IMCLOSE,
%   IMOPEN, STREL, FSPECIAL, RAD2DEG.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Expected radius, area, and position of the reference dot
refRadius = 0.5; % mm
refArea = pi*(refRadius*arena.scale)^2; 
refPosition = [arena.position(1)+0.05*arena.position(3),...
    arena.position(2)+0.95*arena.position(4)]; % bottom-left corner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get threshold for image
thresholdBW = adaptthresh(image,0.6);

% Get the inverse image of the frame
imageInverse = imcomplement(double(image));

% Narrow - broad averaging filters
imageFiltered = imfilter(imageInverse,fspecial('average',3)) - ...
    imfilter(imageInverse,fspecial('average',40),'Replicate');

% Binarize image and mask out arena
imageBW = imbinarize(imageFiltered,thresholdBW*3).*...
    ~(bwdist(arena.mask) <= 10);

% Fill holes and find difference to get just the interior
imageFilled = imfill(imageBW,'holes');
imageDiff = imageFilled - imageBW;

% Get Image Objects
imageLabeled = bwlabel(imageDiff);
imageObjects = regionprops(imageLabeled,'Area','Circularity',...
    'Eccentricity','Centroid');

% Find reference object given circularity, expected area, and position
objectPositions = vertcat(imageObjects.Centroid);
referenceObject = find([imageObjects.Circularity] > 0.4 & ...
    [imageObjects.Area] > refArea*0.12 & ...
    [imageObjects.Area] < refArea*1.5 & ...
    abs(objectPositions(:,1)' - refPosition(1)) < arena.position(3)/10 & ...
        abs(objectPositions(:,2)' - refPosition(2)) < arena.position(3)/10);

% If more than one reference object, use one closest to expected position
if length(referenceObject) > 1
    objectPositions = vertcat(imageObjects(referenceObject).Centroid);
    [~,newReferenceObject] = min(sqrt((objectPositions(:,1)' - refPosition(1)).^2 + ...
        (objectPositions(:,2)' - refPosition(2)).^2));
    if length(newReferenceObject) == 1
        referenceObject = referenceObject(newReferenceObject);
    end
end

% If no reference object found, try morphological closing
if length(referenceObject) ~= 1
    imageClosed = bwdist(~(bwdist(imageBW) <= 3)) >= 3;
    
    % Fill holes and find difference to get just the interior
    imageFilled = imfill(imageClosed,'holes');
    imageDiff = imageFilled - imageClosed;
    
    % Get Image Objects
    imageLabeled = bwlabel(imageDiff);
    imageObjects = regionprops(imageLabeled,'Area','Circularity',...
        'Eccentricity','Centroid');
    
    % Find reference object given circularity, area, position, eccentricity
    objectPositions = vertcat(imageObjects.Centroid);
    referenceObject = find([imageObjects.Circularity] > 0.4 & ...
        [imageObjects.Area] > refArea*0.05 & ...
        [imageObjects.Area] < refArea*1.5 & ...
        abs(objectPositions(:,1)' - refPosition(1)) < arena.position(3)/10 & ...
        abs(objectPositions(:,2)' - refPosition(2)) < arena.position(3)/10 & ...
        [imageObjects.Eccentricity] < 0.7);

    [~,refInd] = max(vertcat(imageObjects(referenceObject).Area));
    referenceObject = referenceObject(refInd);
end

% If no reference object found still, try inverting the image
if length(referenceObject) ~= 1
    imageOpened = bwdist(bwdist(~imageBW) >=2) >= 2;
    
    % Get Image Objects
    imageLabeled = bwlabel(imageOpened);
    imageObjects = regionprops(imageLabeled,'Area','Circularity',...
        'Eccentricity','Centroid');
    
    % Find reference object given circularity, area, position, eccentricity
    objectPositions = vertcat(imageObjects.Centroid);
    referenceObject = find([imageObjects.Circularity] > 0.4 & ...
        [imageObjects.Area] > refArea*0.05 & ...
        [imageObjects.Area] < refArea*1.5 & ...
        abs(objectPositions(:,1)' - refPosition(1)) < arena.position(3)/10 & ...
        abs(objectPositions(:,2)' - refPosition(2)) < arena.position(3)/10 & ...
        [imageObjects.Eccentricity] < 0.7);
end
    
% Get reference dot mask and centroid position
orientation.mask = imageLabeled == referenceObject;
orientation.center = imageObjects(referenceObject).Centroid;

% Calculate the difference from the center of the arena to the center of
% the reference dot
orientation.distance = sqrt((orientation.center(1) - arena.center(1)).^2 + ...
    (orientation.center(2) - arena.center(2)).^2);

% Translate the arena and reference masks
arenaTrans = imtranslate(arena.mask+orientation.mask,arena.offset,'FillValues',0);

% Refind the centroids
imageObjects = regionprops(bwlabel(arenaTrans),'Centroid','Area');
[~,arenaObject] = max([imageObjects.Area]);
[~,refObject] = min([imageObjects.Area]);
centers = vertcat(imageObjects.Centroid);

% Calculate the angle that the pipetting template needs to be rotated by
% (degrees counter-clockwise rotation)
orientation.angle = rad2deg(atan((centers(arenaObject,2) - centers(refObject,2))/...
    (centers(refObject,1) - centers(arenaObject,1)))) - 45;

end