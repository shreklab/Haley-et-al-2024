function [threshold] = getLawnsThreshold(image,arenaMask,lawnRadius)
% [threshold] = GETLAWNSTHRESHOLD(image,arenaMask,lawnRadius)
%
%   GETLAWNSTHRESHOLD creates an estimate of the locations of bacterial 
%   lawns in the behavioral arena given a frame from the video (ideally 
%   without worms) and information about the expected position of the lawns
%   using image processing and morphological operations. GETLAWNSTHRESHOLD
%   expects that the image used is not necessarily from the experimental
%   video itself and enables translation, rotation, and scaling of the lawn
%   mask to match the experimental video.
%
%   INPUTS:
%       - image [1024x1024 uint8]: image of the experiment or lawn video;
%           often the first frame of the contrast video being analyzed
%       - arenaMask [1024x1024 logical]: mask of the arena (1 = inside 
%           arena); if no arena, use false(size(image))
%       - lawnRadius [double]: radii of the expected lawn(s) in pixels 
%           (e.g. 24.5)
%
%   OUTPUTS:
%       - threshold [struct]: a structure containing the following fields
%           - mask [maskSize(1)xmaskSize(2) logical]: mask of the lawn
%               locations (1 = lawn)
%           - labeled [1024x1024 double]: image of the lawn locations with 
%               each lawn uniquely labeled with a number in the set 
%               [1,2,...,numLawn]
%           - image [1024x1024 double]: filtered image
%           - centers [numLawnx2 double]: X,Y positions of each lawn in
%               pixels with row index corresponding to the label(s) in 
%               labeled
%           - radii [numLawnx1 double]: radii of each lawn in pixels with
%               row index corresponding to the label(s) in labeled
%
%   Written 4/13/2022 by Jess Haley in MATLAB R2021a.
%
%   See also ANALYZEFORAGING, GETINVISIBLELAWNS, BWLABEL, REGIONPROPS,
%   GRAYTHRESH, IMBINARIZE, IMCOMPLEMENT, IMFILTER, IMERODE, IMFILL,
%   IMTRANSLATE, IMRESIZE, IMCROP, IMROTATE, PADARRAY, CONV2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get expected lawn area
lawnArea = pi*lawnRadius^2;

% Get the inverse image of the  first frame
firstFrame = double(image);
imageInverse = imcomplement(firstFrame);

% Get threshold for image
thresholdBW = graythresh(firstFrame);

% Narrow - broad averaging filters
imageFiltered = imfilter(imageInverse,fspecial('average',4)) - ...
    imfilter(imageInverse,fspecial('average',60));

% Binarize image and get objects
imageBW = imbinarize(imageFiltered,thresholdBW*3);
if sum(imageBW,'all')/numel(image) < 0.02
    imageBW = imbinarize(imageFiltered,'adaptive','Sensitivity',0.6);
end
imageBW = imerode(imageBW,strel('disk',1)).*arenaMask;
imageLabeled = bwlabel(imageBW);
if max(imageLabeled,[],'all') > 1000
    labeledArea = histcounts(imageLabeled,(0:max(imageLabeled,[],'all'))+0.5);
    largeObjects = find(labeledArea > prctile(labeledArea,100*(1- (1000/max(imageLabeled,[],'all')))));
    imageLabeled = bwlabel(ismember(imageLabeled,largeObjects));
end
imageObjects = regionprops(imageLabeled,'FilledArea','Centroid','Circularity');

% Get objects with similar filled area to expected lawn(s)
lawnObjects = find([imageObjects.FilledArea] > lawnArea*0.7 & ...
    [imageObjects.FilledArea] < lawnArea*3 & ...
    [imageObjects.Circularity] > 0.04);

% Get image objects
lawnMask = false(size(imageBW));
for i = 1:length(lawnObjects)
    lawnMask(imageLabeled == lawnObjects(i)) = true;
end
lawnImage = imfill(lawnMask,'holes');

% Check if enough contrast to try pure thresholding (w/o filtering)
fluorescent = sign(mean(firstFrame(lawnImage),'all') - mean(firstFrame(arenaMask & ~lawnImage),'all'));
valueRange = prctile(firstFrame(lawnImage),[5 99.9]);
if diff(valueRange) > 10 && length(lawnObjects) == 1 && diff(valueRange) < 128
    lastMask = lawnImage; threshMask = lawnImage;
    updateVal = fluorescent;
    while abs(sum(lastMask,'all') - sum(threshMask,'all'))/sum(lawnImage,'all') < 0.2
        updateVal = updateVal - fluorescent*ceil(diff(valueRange)/200);
        lastMask = threshMask;
        threshMask = firstFrame < valueRange(2) + updateVal;
        threshLabeled = bwlabel(threshMask);
        threshMask = threshLabeled == median(threshLabeled(lawnImage));
    end
    lawnImage = imfill(lawnImage | lastMask,'holes');
end

% Blur image via convolution to smooth masks
nublessImage = bwdist(bwdist(~lawnImage)>=10)<=10; % remove nubs
lawnImage = (lawnImage - nublessImage == 0) & (lawnImage > 0);
windowSize = 10;
kernel = ones(windowSize)/windowSize^2;
blurryImage = conv2(single(lawnImage),kernel,'same');
smoothedImage = blurryImage > 0.5;

% Recheck lawn objects
imageLabeled = bwlabel(smoothedImage);
imageObjects = regionprops(imageLabeled,'FilledArea','Centroid','Circularity');
lawnObjects = find([imageObjects.FilledArea] > lawnArea*0.7 & ...
    [imageObjects.FilledArea] < lawnArea*3 & ...
    [imageObjects.Circularity] > 0.5);
lawnMask = false(size(imageBW));
for i = 1:length(lawnObjects)
    lawnMask(imageLabeled == lawnObjects(i)) = true;
end

% Label lawn objects and get centers and radii
threshold = getLawnProperties(lawnMask,imageFiltered);

end

