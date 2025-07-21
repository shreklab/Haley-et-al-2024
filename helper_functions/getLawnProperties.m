function [imageStruct] = getLawnProperties(mask,image)
% [imageStruct] = GETLAWNPROPERTIES(mask,image)
%
%   GETLAWNPROPERTIES takes a binary mask of regions of interest (e.g.,
%   bacterial lawns), labels each distinct region, and computes their
%   geometric properties.
%
%   INPUTS:
%       - mask [MxN logical]: A binary mask where true values indicate the
%           locations of the regions to be analyzed.
%       - image [MxN num]: (Optional) The original image associated with
%           the mask.
%
%   OUTPUTS:
%       - imageStruct [struct]: A structure containing the following fields:
%           - image [MxN num]: The original image, if provided.
%           - mask [MxN logical]: The input binary mask of all regions.
%           - labeled [MxN double]: An image where each distinct region
%               from the mask is uniquely labeled with an integer.
%           - centers [numLawnx2 double]: The [X,Y] pixel coordinates for
%               the centroid of each labeled region.
%           - radii [numLawnx1 double]: The radius of each labeled region in
%               pixels, calculated from its area.
%           - circularity [numLawnx1 double]: The circularity of each
%               labeled region, where 1.0 indicates a perfect circle.
%
%   Written 3/26/2024 by Jess Haley in MATLAB R2024a.
%
%   See also GETLAWNSTHRESHOLD, GETLAWNSTEMPLATE, GETLAWNSFILTER, 
%   GETLAWNSCIRCLE, BWLABEL, REGIONPROPS.

if nargin == 2
    imageStruct.image = image;
end
imageStruct.mask = mask;
imageStruct.labeled = uint8(bwlabel(imageStruct.mask));
imageObjects = regionprops(imageStruct.labeled,'Centroid','Area','Circularity');
imageStruct.centers = vertcat(imageObjects.Centroid);
imageStruct.radii = sqrt(vertcat(imageObjects.Area)/pi);
imageStruct.circularity = vertcat(imageObjects.Circularity);

end