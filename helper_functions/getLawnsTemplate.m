function [template] = getLawnsTemplate(maskSize,scale,arenaDiameter,lawnDiameter,lawnSpacing)
% [template] = GETLAWNSTEMPLATE(maskSize,scale,arenaDiameter,lawnDiameter,lawnSpacing)
%
%   GETLAWNSTEMPLATE creates an estimate of the locations of bacterial 
%   lawns in the behavioral arena given information about size and spacing 
%   of a single lawn or an isometric grid of lawns.
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
%       - template [struct]: a structure containing the following fields
%           - mask [maskSize(1)xmaskSize(2) logical]: mask of the lawn
%               locations (1 = lawn)
%           - labeled [maskSize(1)xmaskSize(2) double]: image of the lawn
%               locations with each lawn uniquely labeled with a number in
%               the set [1,2,...,numLawn]
%           - centers [numLawnx2 double]: X,Y positions of each lawn in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - radii [numLawnx1 double]: radii of each lawn in pixels with
%               row index corresponding to the label(s) in labeled
%
%   Written 4/13/2022 by Jess Haley in MATLAB R2021a.
%
%   See also ANALYZEFORAGING, GETINVISIBLELAWNS, BWLABEL, REGIONPROPS,
%   REPMAT, FLOOR, MOD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If only one lawn, center on the mask
if lawnSpacing == 0
    % Get center and radius
    centers = maskSize/2;
    radii = lawnDiameter*scale/2;

% If multiple lawns, create an isometric grid
else
    
    % Get the center-to-center distance between lawns in X and Y
    lawnDistY = lawnSpacing; % center-to-center distance
    lawnDistX = lawnDistY*sqrt(3)/2; % equilateral triangle
    
    % Calculate the number of lawns that can fit in the arena
    lawnRepeats = floor(2*arenaDiameter/lawnDistY); % # of lawns that can fit
    lawnRepeats = lawnRepeats + mod(lawnRepeats,2) - 1; % force to be odd
    
    % Get the X and Y positions of the lawn centers
    lawnPosX = lawnDistX.*[1:lawnRepeats];
    lawnPosX = repmat(lawnPosX,lawnRepeats,1);
    lawnPosY = lawnDistY/2.*[1:lawnRepeats];
    lawnPosY = repmat(lawnPosY,lawnRepeats,1) + 2*lawnPosY';
    
    % Center the grid on the origin [0,0]
    origin = [mean(lawnPosX,'all'),mean(lawnPosY,'all')];
    lawnPosX = lawnPosX - origin(1);
    lawnPosY = lawnPosY - origin(2);
    
    % Find lawns within arena boundary
    bounds = arenaDiameter(1)/2 - lawnSpacing/2;
    in = sqrt(lawnPosX.^2 + lawnPosY.^2) <= bounds;
    lawnPosX = lawnPosX(in); lawnPosY = lawnPosY(in); % mm
    numLawns = length(lawnPosX);
    
    % Convert to pixels (and swap x and y)
    xPix = lawnPosY*scale + maskSize(2)/2;
    yPix = lawnPosX*scale + maskSize(1)/2; % pix
    
    % Get centers and radii
    centers = [xPix,yPix];
    radii = repmat(scale*lawnDiameter/2,numLawns,1);
end

% Create lawn mask from centers and radii
X = repmat(1:maskSize(1),maskSize(2),1)';
Y = repmat(1:maskSize(2),maskSize(1),1);
mask = false(maskSize);
for numROI = 1:length(radii)
    distances = sqrt((X-centers(numROI,2)).^2 + ...
        (Y-centers(numROI,1)).^2);
    mask(distances <= radii(numROI)) = 1;
end

% Label lawn objects and get centers and radii
template = getLawnProperties(mask);

end

