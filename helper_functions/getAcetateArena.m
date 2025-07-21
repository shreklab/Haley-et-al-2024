function [arena] = getAcetateArena(image,arenaSize)
% [arena] = GETACETATEARENA(image,arenaSize,flag)
%
%   GETACETATEARENA takes in an image containing a piece of acetate with a 
%   large hole cut out of it (i.e. an arena) to contain worm behavior and 
%   returns a mask of the arena and its properties using image processing 
%   and morphological operations.
%
%   INPUTS:
%       - image [1024x1024 uint8]: image of the acetate arena; often the
%           first frame of the video being analyzed
%       - arenaSize [double OR 1x2 double]: number giving the diameter or
%           length and width of the arena in mm (e.g. 30 OR [16,25])
%
%   OUTPUTS:
%       - arena [struct]: a structure containing the following fields
%           - mask [1024x1024 logical]: mask of the arena (1 = inside arena)
%           - outline [1024x1024 logical]: mask of the arena perimeter 
%               (1 = perimeter of arena only)
%           - center [1x2 double]: X,Y location of the arena centroid in
%               pixels
%           - offset [1x2 double]: X,Y offset of the arena centroid from
%               the image center in pixels
%           - area [double]: area of the arena in pixels^2
%           - position [1x4 double]: bounding box of the arena [X,Y,W,H] in
%               pixels
%           - circularity [double]: circularity of the arena (e.g. 1 =
%               perfect circle)
%           - majorAxisLength [double]: length of the longest axis of the
%               arena ellipse in pixels
%           - minorAxisLength [double]: length of the shortest axis of the
%               arena ellipse in pixels
%           - meanIntensity [double]: mean value of all the pixels in
%               image within the arena mask
%           - minIntensity [double]: minimum value of all the pixels in
%               image within the arena mask
%           - maxIntensity [double]: maximum value of all the pixels in
%               image within the arena mask
%           - shape [double]: confirmation that the arena mask was found
%               correctly; either 'circle' or 'error' based on circularity
%           - diameter [double]: arena diameter given in mm
%           - scale [double]: calculated scale of the arena given the 
%               measured area and given diameter in pixels/mm
%    
%   Written 4/13/2022 by Jess Haley in MATLAB R2021a.
%
%   See also ANALYZEFORAGING, GETACETATEORIENTATION, REGIONPROPS, BWLABEL, 
%   BWMORPH, ADAPTTHRESH, IMBINARIZE, IMDILATE, IMFILL, STREL, CONV2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if arena is circular or rectangular
if length(arenaSize) == 1
    arenaShape = 'circle';
    arenaDiameter = arenaSize;
else
    arenaShape = 'rectangle';
end

% Get threshold for image
thresholdBW = adaptthresh(image,0.6);

% Binarize image
imageBW = imbinarize(image,thresholdBW);
imageBW = imageBW + bwdist(imageBW) == 1;

% Get largest image object 
imageLabeled = bwlabel(imageBW);
imageObjects = regionprops(imageLabeled,'FilledArea','Circularity');
if strcmp(arenaShape,'circle')
    arenaObject = find([imageObjects.FilledArea] > 0.35*numel(image) & ...
        [imageObjects.Circularity] > 0.55);
    if imageObjects(arenaObject).Circularity < 0.99 && arenaSize == 9
        imageBW = imbinarize(image);s
        imageBW = imageBW + bwdist(imageBW) == 1;
        imageLabeled = bwlabel(imageBW);
        imageObjects = regionprops(imageLabeled,'FilledArea','Circularity');
        arenaObject = find([imageObjects.FilledArea] > 0.35*numel(image) & ...
            [imageObjects.Circularity] > 0.6);
    end
else
    arenaObject = find([imageObjects.FilledArea] > 0.35*numel(image) & ...
        [imageObjects.FilledArea] < 0.75*numel(image));
end

% If no arena object detected, try lowering the threshold
if isempty(arenaObject)
    thresholdBW = adaptthresh(image,0.59);
    imageBW = imbinarize(image,thresholdBW);
    imageBW = imageBW + bwdist(imageBW) == 1;
    imageLabeled = bwlabel(imageBW);
    imageObjects = regionprops(imageLabeled,'FilledArea','Circularity');
    arenaObject = find([imageObjects.FilledArea] > 0.5*numel(image) & ...
        [imageObjects.Circularity] > 0.2);
end

% Fill arena mask
imageFilled = imfill(imageLabeled == arenaObject,'holes');
imageFilled = imfill(imageFilled + bwdist(imageFilled) == 1,'holes');

% Smooth arena mask
windowSize = 10;
kernel = ones(windowSize)/windowSize^2;
blurryImage = conv2(single(imageFilled),kernel,'same');
arena.mask = imfill(blurryImage > 0.5,'holes');

% Get arena properties
arenaProps = regionprops(arena.mask,image,'Centroid','Area','Circularity',...
    'Eccentricity','MajorAxisLength','MinorAxisLength','WeightedCentroid',...
    'MeanIntensity','MaxIntensity','MinIntensity','BoundingBox','Orientation');
arena.outline = bwmorph(arena.mask,'remove');
arena.center = flip(arenaProps.Centroid);
arena.offset = arena.center - size(arena.mask)/2;
arena.area = arenaProps.Area;
arena.position = arenaProps.BoundingBox;
arena.circularity = arenaProps.Circularity;
arena.majorAxisLength = arenaProps.MajorAxisLength;
arena.minorAxisLength = arenaProps.MinorAxisLength;
arena.meanIntensity = arenaProps.MeanIntensity;
arena.minIntensity = arenaProps.MinIntensity;
arena.maxIntensity = arenaProps.MaxIntensity;
if arena.circularity > 0.9
    arena.shape = 'circle';
    arena.diameter = arenaDiameter; % mm
    arena.scale = sqrt(arena.area/(pi*(0.5*arenaDiameter)^2)); % pixels/mm
elseif  arena.area/prod(arena.position([3,4])) > 0.89
    arena.shape = 'rectangle';
    arena.size = arenaSize; % mm
    arena.scale = sqrt(arena.area/prod(arenaSize)); % pixels/mm
    arena.angle = arenaProps.Orientation;
else
    arena.shape = 'error';
end

end