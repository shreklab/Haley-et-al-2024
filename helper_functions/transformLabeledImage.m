function [transformedStruct] = transformLabeledImage(imageStruct,translate,rotate,resize,properties)
% [template] = TRANSFORMLABELEDIMAGE(labeledImage,translate,rotate,resize,properties)
%
%   TRANSFORMLABELEDIMAGE takes in an image and it's properties as a
%   imageStruct (created by functions such as CREATEHEXAGONALGRID and
%   CREATEDIAMONDGRID) and transforms the image with given translation,
%   rotation, and resizing. The transformed image is then returned as new
%   structure with corresponding properties.
%
%   INPUTS:
%       - imageStruct [struct]: a structure containing the following 
%               fields (labeled is required; the rest are optional):
%           - labeled [MxN double]: image of the ROI locations uniquely 
%               labeled with a number in the set [1,2,...,numROI]
%           - mask [MxN logical]: mask of all ROI locations (1 = on)
%           - centers [numROIx2 double]: X,Y positions of each ROI in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - radii [numROIx1 double]: radii of each ROI in pixels with
%               row index corresponding to the label(s) in labeled
%           - area [numROIx1 double]: area of each ROI in pixels with
%               row index corresponding to the label(s) in labeled
%       - translate [1x2 double]: the number of pixels in X,Y that the 
%           image needs to be translated by (e.g. [20.7, -4.2])
%       - rotate [double]: the number of degrees the image needs to be 
%           rotated by (e.g. -3.1)
%       - resize [double]: the scale at which the image needs to be resized
%           by (e.g. 1.2)
%       - properties [cell array]: list of properties to be returned (i.e.
%           {'centers','area','radius'};
%
%   OUTPUTS:
%       - transformedStruct [struct]: a structure containing the following 
%               fields (labeled is required; the rest are optional):
%           - labeled [MxN double]: image of the ROI locations uniquely 
%               labeled with a number in the set [1,2,...,numROI]
%           - mask [MxN logical]: mask of all ROI locations (1 = on)
%           - centers [numROIx2 double]: X,Y positions of each ROI in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - radii [numROIx1 double]: radii of each ROI in pixels with
%               row index corresponding to the label(s) in labeled
%           - area [numROIx1 double]: area of each ROI in pixels with
%               row index corresponding to the label(s) in labeled
%
%   Written 5/2/2022 by Jess Haley in MATLAB R2022a.
%
%   See also ANALYZEFORAGING, CREATEHEXAGONALGRID, CREATEDIAMONDGRID,
%   IMTRANSLATE, IMRESIZE, IMCROP, PADARRAY, IMROTATE, REGIONPROPS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get image and fields
image = imageStruct.labeled;
properties = fields(imageStruct);

% Translate, resize, and rotate image
translatedImage = imtranslate(image,translate,'nearest','FillValues',0);
if resize ~= 1
    resizedImage = imresize(translatedImage,resize,'nearest');
    if size(resizedImage,1) > size(translatedImage)
        croppedImage = imcrop(resizedImage,...
            centerCropWindow2d(size(resizedImage),size(translatedImage)));
    else
        croppedImage = padarray(resizedImage,...
            (size(translatedImage)-size(resizedImage))/2,'both');
    end
else
    croppedImage = translatedImage;
end

labeled = imrotate(croppedImage,rotate,'nearest','crop');

% If mask is a property of the original struct, find it
if sum(strcmp(properties,'mask')) > 0
        transformedStruct.mask = labeled > 0;
end

transformedStruct.labeled = labeled;
imageObjects = regionprops(transformedStruct.labeled,'Centroid','Area','Circularity');

% If centers are a property of the original struct, find them
if sum(strcmp(properties,'centers')) > 0
        transformedStruct.centers = vertcat(imageObjects.Centroid);
end

% If radii are a property of the original struct, find them
if sum(strcmp(properties,'radii')) > 0
        transformedStruct.radii = sqrt(vertcat(imageObjects.Area)/pi);
end

% If area is a property of the original struct, find it
if sum(strcmp(properties,'area')) > 0
        transformedStruct.area = vertcat(imageObjects.Area);
end

if sum(strcmp(properties,'circularity')) > 0
        transformedStruct.circularity = vertcat(imageObjects.Circularity);
end

end