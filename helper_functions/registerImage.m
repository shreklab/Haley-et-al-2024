function [imageC,tform] = registerImage(imageA,imageB,pointsA,pointsB)
%[imageC,tform] = REGISTERIMAGE(imageA,imageB,pointsA,pointsB)
%   
%   REGISTERIMAGE takes in two images, detects SURF features of each, and 
%   performs a rigid transformation of imageB to imageA. The transformed 
%   imageB is then returned along with the corresponding 2D geometric 
%   transformation.
%
%   INPUTS
%       - imageA [MxN array]: reference image (e.g., first frame of the
%           behavior video)
%       - imageB [MxN array]: image needing transformation (e.g., first
%           frame of the lawn contrast video)
%       - pointsA [SURFPoints]: object containing an array of SURF point
%           coordinates for imageA (optional)
%       - pointsB [SURFPoints]: object containing an array of SURF point
%           coordinates for imageB (optional)
%
%   OUTPUTS
%       - imageC [MxN array]: imageB after transformation
%       - tform [rigidtform2d]: object containing the 2D transformation of
%           imageB into imageC
%
%   Written 2/24/2024 by Jess Haley in MATLAB R2023b.
%
%   See also GETINVISIBLELAWNS, DETECTSURFFEATURES, EXTRACTFEATURES,
%   MATCHFEATURES, ESTGEOTFORM2D, IMWARP, SHOWMATCHEDFEATURES, IMSHOWPAIR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    % Detect features
    pointsA = []; pointsB = [];
    threshA = 55; threshB = 50;
    while size(pointsA,1) < 400 && threshA >= 10
        pointsA = detectSURFFeatures(imageA,'MetricThreshold',threshA);
        threshA = threshA - 10;
    end
    while size(pointsB,1) < 1500 && threshB >= 5
        pointsB = detectSURFFeatures(imageB,'MetricThreshold',threshB);
        threshB = threshB - 5;
    end
    [featA,validA] = extractFeatures(imageA,pointsA);
    [featB,validB] = extractFeatures(imageB,pointsB);
    [indexPairs] = matchFeatures(featA,featB,'MatchThreshold',100,'MaxRatio',1);
    % if size(indexPairs,1) <2
    %     pointsA = detectMSERFeatures(imageA,'RegionAreaRange',[30 round(numel(imageA)*0.90)]);
    %     pointsB = detectMSERFeatures(imageB,'RegionAreaRange',[30 round(numel(imageA)*0.90)]);
    %     [featA,validA] = extractFeatures(imageA,pointsA);
    %     [featB,validB] = extractFeatures(imageB,pointsB);
    %     indexPairs = matchFeatures(featA,featB,'MatchThreshold',100,'MaxRatio',1);
    % end
    matchedA = validA(indexPairs(:,1));
    matchedB = validB(indexPairs(:,2));
    [tform,inlierIdx] = estgeotform2d(matchedB,matchedA,'similarity','MaxNumTrials',100000);
    % if abs(tform.Scale - 1) > 0.2
    %     pointsA = detectMSERFeatures(imageA,'RegionAreaRange',[30 round(numel(imageA)*0.90)]);
    %     pointsB = detectMSERFeatures(imageB,'RegionAreaRange',[30 round(numel(imageA)*0.90)]);
    %     [featA,validA] = extractFeatures(imageA,pointsA);
    %     [featB,validB] = extractFeatures(imageB,pointsB);
    %     indexPairs = matchFeatures(featA,featB,'MatchThreshold',100,'MaxRatio',1);
    %     matchedA = validA(indexPairs(:,1));
    %     matchedB = validB(indexPairs(:,2));
    %     [tform,inlierIdx] = estgeotform2d(matchedB,matchedA,'similarity');
    % end
    inlierA = matchedA(inlierIdx);
    inlierB = matchedB(inlierIdx);
elseif nargin == 4
    distThreshold = min(pdist(pointsB))/2;
    distPoints = pdist2(pointsA,pointsB);
    [indA,indB] = find(distPoints < distThreshold);
    matchedA = pointsA(indA,:);
    matchedB = pointsB(indB,:);
    [tform,inlierIdx] = estgeotform2d(matchedB,matchedA,'rigid',...
        'MaxDistance',distThreshold);
    inlierA = matchedA(inlierIdx,:);
    inlierB = matchedB(inlierIdx,:);
end

imageC = imwarp(imageB,tform,OutputView=imref2d(size(imageA)));

figure
subplot(121); showMatchedFeatures(imageA,imageB,inlierA,inlierB)
subplot(122); imshowpair(imageA,imageC)

end