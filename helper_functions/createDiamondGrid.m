function [diamond,triangle] = createDiamondGrid(maskSize,scale,arenaDiameter,lawnDiameter,lawnSpacing)
% [diamond,triangle] = CREATEDIAMONDGRID(maskSize,scale,arenaDiameter,lawnDiameter,lawnSpacing)
%
%   CREATEDIAMONDGRID creates a grid of diamonds and corresponding 
%   triangles (2 per diamond) given information about their size and 
%   spacing.
%
%   INPUTS:
%       - maskSize [1x2 double]: the width and height of the desired output
%           image (e.g. [1024, 1024])
%       - scale [double]: the image scaling in pixels per mm (e.g. 32)
%       - arenaDiameter [double]: the diameter of the arena in mm (e.g. 30)
%       - lawnDiameter [double]: the diameter of the lawn(s) in mm (e.g.
%           1.5)
%       - lawnSpacing [double]: the center-to-center distance between lawns
%           in mm (e.g. 6); use 0 for a single lawn
%
%   OUTPUTS:
%       - diamond [struct]: a structure containing the following fields
%           - labeled [maskSize(1)xmaskSize(2) double]: image of the grid
%               locations uniquely labeled with a number in the set 
%               [1,2,...,numGrid]
%           - centers [numGridx2 double]: X,Y positions of each grid in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - area [numGridx1 double]: area of each grid in pixels with
%               row index corresponding to the label(s) in labeled
%       - triangle [struct]: a structure containing the following fields
%           - labeled [maskSize(1)xmaskSize(2) double]: image of the grid
%               locations uniquely labeled with a number in the set 
%               [1,2,...,numGrid]
%           - centers [numGridx2 double]: X,Y positions of each grid in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - area [numGridx1 double]: area of each grid in pixels with
%               row index corresponding to the label(s) in labeled
%
%   Written 4/13/2022 by Jess Haley in MATLAB R2022a.
%
%   See also ANALYZEFORAGING, GETLAWNSTEMPLATE, IMERODE, STREL, BWLABEL,
%   REGIONPROPS, IMDILATE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get lawn template (isometric grid)
template = getLawnsTemplate(maskSize,scale,2*arenaDiameter,lawnDiameter,lawnSpacing/2);

% Define slope and x values 
m = sqrt(3)/3;
x = (-2*arenaDiameter:lawnSpacing:2*arenaDiameter) + lawnSpacing/2;

% Create grids of X and Y locations
imageSize = maskSize/scale; % mm [width x height]
X = repmat(linspace(-imageSize(1)/2,imageSize(1)/2,maskSize(1)),maskSize(2),1)';
Y = repmat(linspace(-imageSize(2)/2,imageSize(2)/2,maskSize(2)),maskSize(1),1);

% Create mask of equilateral triangles anti-aligned with template
diamondMask = zeros(maskSize);
for i = 1:length(x)
    mask1 = zeros(maskSize);
    mask1(X./(Y-x(i)) >= m) = 1;
    mask1(X./(Y-x(i)) <= -m) = 2;
    mask2 = zeros(maskSize);
    mask2(Y-lawnSpacing/2 >= x(i)) = 1;
    diamondMask = diamondMask + mask1 + mask2;
end

% Identify individual triangles and label uniquely
values = unique(diamondMask(:));
gridLabeled = zeros(maskSize);
gridCount = 0;
for i = 1:length(values)
    imageBW = imerode(diamondMask == values(i),strel('square',2));
    imageLabeled = bwlabel(imageBW);
    numObjects = length(unique(imageLabeled(:)));
    for j = 1:numObjects
        gridLabeled(imageLabeled == j) = j + gridCount;
    end
    gridCount = gridCount + numObjects - 1;
end

% Remove partial triangles and triangles not in the arena
gridObjects = regionprops(gridLabeled,'Area');
triangleObjects = find([gridObjects.Area] > 0.9*median([gridObjects.Area]));
% gridObjects = regionprops(gridLabeled,arena.mask,'Area','MeanIntensity');
% triangleObjects = find([gridObjects.Area] > 0.9*median([gridObjects.Area]) & ...
%     [gridObjects.MeanIntensity] > 0);
trianglesLabeled = zeros(maskSize);
for i = 1:length(triangleObjects)
    trianglesLabeled(gridLabeled == triangleObjects(i)) = i;
end

% Group triangles by proximity to lawn centers
gridObjects = regionprops(trianglesLabeled,'Centroid','Area');
gridCenters = vertcat(gridObjects.Centroid);
templateCenters = template.centers;
diamond.labeled = zeros(maskSize);
triangle.labeled = zeros(maskSize);
for i = 1:size(templateCenters,1)
    dist = sqrt((gridCenters(:,1)-templateCenters(i,1)).^2+...
        (gridCenters(:,2)-templateCenters(i,2)).^2);
    ind = find(dist < 0.2*lawnSpacing*scale);
    for j = 1:length(ind)
        mask = imdilate(trianglesLabeled == ind(j),strel('diamond',2));
        diamond.labeled(mask) = double(length(ind) == 2) * ...
            (max(diamond.labeled,[],'all') + double(j == 1));
        triangle.labeled(mask) = max(triangle.labeled,[],'all') + 1;
    end
end

% Label grid objects and get centers and area
diamondObjects = regionprops(diamond.labeled,'Centroid','Area');
diamond.centers = vertcat(diamondObjects.Centroid);
diamond.area = vertcat(diamondObjects.Area);
triangleObjects = regionprops(triangle.labeled,'Centroid','Area');
triangle.centers = vertcat(triangleObjects.Centroid);
triangle.area = vertcat(triangleObjects.Area);

end

