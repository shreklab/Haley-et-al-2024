function [data] = analyzeWormLabTracks(info,plateInd,folder,bodyPart)
%[data] = ANALYZEWORMLABTRACKS(info,expToAnalyze)
%
%   ANALYZEWORMLABTRACKS takes mid-point data exported from WormLab and 
%   computes numerous metrics related to the worm's position, velocity,
%   and turning behavior as well as it's location relative to the arena and
%   lawns.
%
%   INPUTS:
%       - info [table]: contains metadata about the experiments; see
%           GETEXPERIMENTINFO for more information about the contents of 
%           this variable
%       - plateInd [double OR logical array]: plate id to analyze
%       - folder [char]: directory to load the files from
%       - bodyPart [str]: specify the body-part to be analyzed (if
%           applicable)
%   OUTPUTS:
%       - data [table]: a table containing the following columns
%           - expNum [double]: unique numerical identifier for an
%               experiment
%           - plateNum [double]: unique numerical identifier for a single
%               assay plate
%           - videoNum [double]: video file number (applicable if multiple
%               files were saved for that plate)
%           - wormNum [double]: unique numerical identifier for an
%               individual worm
%           - frameNum [double]: frame number
%           - time [double]: time stamp of the frame in seconds
%           - xPosition [double]: X position of the worm's midpoint in 
%               pixels at the given frame (e.g. 611.2309); if NaN, worm was
%               not located at this position
%           - yPosition [double]: Y position of the worm's midpoint in 
%               pixels at the given frame (e.g. 442.8182); if NaN, worm was
%               not located at this position
%           - bodyPart [str]: body part tracked (i.e. 'head','midpoint', or
%               'tail')
%           - distance [double] distance traveled in one frame in pixels
%           - velocitySmooth [double]: smoothed velocity in microns/second
%               (movmean over a 4 (+/- 2) second time window)
%           - velocityEuclidean [double] Euclidean velocity in
%               microns/second (over a 4 (+/- 2) second time window)
%           - pathAngle [double]: signed path angle in degrees (-180:180) 
%               for a +/- 4 second time window; negative = right turn, 
%               positive = left turn; the larger the magnitude, the greater
%               the reorientation (e.g. -10.2334)
%           - turn [logical]: reorientation events with reorientations of
%               path angle >=150 degrees in a +/- 4 second time window
%               (1 = turn)
%           - distanceNeighbor1 [double]: distance to nearest neighbor in
%               pixels
%           - distanceNeighbor2 [double]: distance to next nearest neighbor
%               in pixels
%           - distanceNeighbor3 [double]: distance to thirds nearest 
%               neighbor in pixels
%           - closestLawnID [double]: gives the numerical identifier of the 
%               lawn an animal is nearest to (e.g. 1)
%           - closestLawnOD600 [double]: gives the OD600 of the lawn an 
%               animal is nearest to (e.g. 10)
%           - distanceLawnEdge [double]: distance from the animal's
%               mid-point to the closest lawn edge in pixels (negative =
%               off lawn, positive = on lawn)
%           - distanceArenaEdge [double]: distance from the animal's
%               mid-point to the arena edge in pixels
%           - nearArena [double]: indicates if an animal is within 1 mm of 
%               the arena edge (1 = near arena)
%           - outOfBounds [double]: indicates that an animal's track was
%               lost at the arena edge and that it is likely hiding at the
%               arena edge or escaped under or on top of the acetate arena
%
%   Written 3/25/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEFORAGING, GETEXPERIMENTINFO, READMATRIX, MOVMEAN,
%   MOVSUM, ATAN2D, FINDPEAKS, PDIST, SQUAREFORM, SORT, REPMAT, RESHAPE,
%   SUB2IND, IMDILATE, IMERODE, STREL, BWDIST.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load Worm Tracks as Exported from WormLab

% Get file name, image size, and scale for this experiment
fileName = info.wormLabFileName{plateInd};
pixels = info.pixels{plateInd};
scale = info.scale(plateInd);

% Load the tracks for this plate
theseTracks = readmatrix(fullfile(folder,fileName),'Range',[6,1]);
if contains(folder,'centerpoints')
    numSections = (size(theseTracks,2)-2)/2;
    if strcmp(bodyPart,'head')
        theseTracks = theseTracks(:,1:4);
    elseif strcmp(bodyPart,'midpoint')
        theseTracks = theseTracks(:,[1:2,2+numSections,3+numSections]);
    elseif strcmp(bodyPart,'tail')
        theseTracks = theseTracks(:,[1:2,1+numSections*2,2+numSections*2]);
    end
end
maxFrame = size(theseTracks,1);
numTracks = (size(theseTracks,2) - 2)/2;

% Add NaN values to the end of the tracks if tracks are shorter
% than the longest tracks
theseTracks = [theseTracks;nan(maxFrame - size(theseTracks,1),size(theseTracks,2))];

% Get frame rate
time = theseTracks(:,2);
dTime = [diff(time(1:2));diff(time)];
frameRate = 1/mean(dTime,'omitnan');

% Store all tracks into a matrix
trackX = nan(maxFrame,numTracks); % all X positions
trackY = nan(maxFrame,numTracks); % all Y positions
trackMissing = false(maxFrame,numTracks);

for k = 1:numTracks

    % Get the X,Y position of the current track
    thisTrack = theseTracks(:,2*k + [1,2]);

    % If missing values are less than 10 seconds long, connect track
    % via linear interpolation
    [thisTrack,noTrack] = fillmissing(thisTrack,'linear','SamplePoints',time,'MaxGap',10);
    trackMissing(:,k) = any(noTrack,2) | thisTrack(:,1) > pixels(2) |...
        thisTrack(:,2) > pixels(1) | any(thisTrack <= 0.5,2) | any(isnan(thisTrack),2);

    % Delete unused columns and pixel values outside of image bounds
    thisTrack(thisTrack(:,1) > pixels(1) |...
        thisTrack(:,2) > pixels(2),:) = NaN;
    thisTrack(thisTrack <= 0.5) = NaN;

    % Add current track to track positions
    trackX(:,k) = thisTrack(:,1);
    trackY(:,k) = thisTrack(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Compute Distance Traveled and Velocity

% Calculate distance traveled at each time step
distance = [nan(1,numTracks);sqrt(diff(trackX).^2 + diff(trackY).^2)];

% Define a smoothing window for calculating distance traveled
winTime = 4; % (+/- winTime/2) second moving window
winSize = round(winTime*frameRate); % frames

% Find and replace outliers
distanceOutliers = (1000/scale)*distance./dTime > 600; % 600 um/s

% Use 4 (+/-2) second moving time window to smooth distance traveled
distanceSmooth = distance; distanceSmooth(distanceOutliers) = NaN;
distanceSmooth = movmean(distanceSmooth,winSize,'omitnan');

% Use 4 (+/-2) second moving time window to calculate Euclidean distance
dX = trackX(1+winSize:end,:) - trackX(1:end-winSize,:);
dY = trackY(1+winSize:end,:) - trackY(1:end-winSize,:);
distanceEuclidean = [nan(ceil(winSize/2),numTracks);...
    sqrt(dX.^2 + dY.^2);nan(floor(winSize/2),numTracks)]./winSize;

% Get velocity (microns/second)
velocitySmooth = (1000/scale)*distanceSmooth./dTime;
velocityEuclidean = (1000/scale)*distanceEuclidean./dTime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Compute Path Angle and Detect Turn Events

% Define a smoothing window for calculating path angle
winTime = 8; %(+/- winTime/2) second moving window
winSize = ceil(winTime*frameRate/2); % frames

% Calculate vector components (winSize long)
dX = [nan(winSize,numTracks);...
    trackX(1+winSize:end,:) - trackX(1:end-winSize,:)];
dY = [nan(winSize,numTracks);...
    trackY(1+winSize:end,:) - trackY(1:end-winSize,:)];

% Start parellel pool
% if isempty(gcp('nocreate'))
%     pool = parpool(feature('numcores'));
% end

% Calculate directional path angle and distance traveled in +/- 4 second
% moving window
% EXPLANATION: https://www.mathworks.com/matlabcentral/answers/180131-how-can-i-find-the-angle-between-two-vectors-including-directional-information
pathAngle = nan(maxFrame,numTracks);
for j = 1:numTracks
    % parfor k = 1:maxFrame-winSize
    for k = 1:maxFrame-winSize
        u = [dX(k,j),dY(k,j)];
        v = [dX(k+winSize,j),dY(k+winSize,j)];
        pathAngle(k,j) = atan2d(u(1)*v(2)-u(2)*v(1),u(1)*v(1)+u(2)*v(2));
    end
end

% Find turning events (> 150 deg)
turn = false(maxFrame,numTracks);
for j = 1:numTracks
    [~,turnLoc] = findpeaks(abs(pathAngle(:,j)),'MinPeakHeight',150,...
        'MinPeakDistance',winSize*2);
    turn(turnLoc,j) = 1;
end
    
% If inadequate speed (less than 30 microns/sec on average), 
% discount turn
% pathLength = sqrt(dX.^2 + dY.^2);
% pathSum = 1000*movmean(pathLength,[winSize winSize],'omitnan')./scale;
% turn(pathSum <= 30) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Compute Worm-to-Worm Distance

% Get distance to 3 nearest neighbors
neighbor = nan(3,maxFrame,numTracks);
for j = 1:maxFrame

    % Interpolate b/w track points before computing neighboring distances

    % Compute pair-wise distance between X,Y positions
    distMatrix = squareform(pdist([trackX(j,:);trackY(j,:)]'));
    distMatrix = sort(distMatrix,2)';
    k = 1;
    while k < numTracks
        neighbor(k,j,:) = distMatrix(k+1,:);
        k = k + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Organize Data into a Table

data = table();

% Add identifiers
data.expNum = repmat(info.expNum(plateInd),maxFrame*numTracks,1);
data.plateNum(:) = info.plateNum(plateInd);
data.videoNum(:) = info.videoNum(plateInd);
[frameNum,trackNum] = find(trackX);
data.wormNum = reshape(info.wormNum{plateInd}(trackNum),[],1);
data.frameNum = frameNum;
data.time = repmat(theseTracks(:,2),numTracks,1);

% Add X,Y position and body part
data.xPosition = reshape(trackX,maxFrame*numTracks,1);
data.yPosition = reshape(trackY,maxFrame*numTracks,1);
data.bodyPart(:) = {bodyPart};

% Add indices where tracks did not exist in WormLab
data.noTrack = reshape(trackMissing,maxFrame*numTracks,1);

% Add distance traveled
data.distance = reshape(distance,maxFrame*numTracks,1);

% Add velocity
data.velocitySmooth = reshape(velocitySmooth,maxFrame*numTracks,1);
data.velocityEuclidean = reshape(velocityEuclidean,maxFrame*numTracks,1);

% Add path angle and turning events
data.pathAngle = reshape(pathAngle,maxFrame*numTracks,1);
data.turn = reshape(turn,maxFrame*numTracks,1);

% Calculate cumulative turns and turning rate with 5 min moving window
% winTime = 5; % 5 min moving window
% winSize = round(60*winTime*frameRate); % frames
% data.turnCum = reshape(cumsum(turn),maxFrame*numTracks,1);
% data.turnRate = reshape(60*frameRate.*movmean(turn,winSize,'omitnan'),...
%     maxFrame*numTracks,1);

% Add distance to 3 nearest neighbors
data.distanceNeighbor1 = reshape(neighbor(1,:,:),maxFrame*numTracks,1);
data.distanceNeighbor2 = reshape(neighbor(2,:,:),maxFrame*numTracks,1);
data.distanceNeighbor3 = reshape(neighbor(3,:,:),maxFrame*numTracks,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Compute Relative Distance to Arena and Lawn(s)

% Convert X,Y positions to index notation
x = round(data.xPosition); x(data.noTrack) = 1;
y = round(data.yPosition); y(data.noTrack) = 1;
indices = sub2ind(pixels,x,y);

% Get distance to nearest lawn edge
mask = flip(info.lawnMask{plateInd})';
mask = -bwdist(mask) + bwdist(~imerode(mask,strel('disk',1)));
data.distanceLawnEdge = double(mask(indices));
data.distanceLawnEdge(data.noTrack) = NaN;

% Get distance to nearest arena edge
mask = flip(info.arenaMask{plateInd})';
mask = bwdist(~mask) - bwdist(~imerode(~mask,strel('disk',1)));
data.distanceArenaEdge = double(mask(indices));
data.distanceArenaEdge(data.noTrack) = NaN;

% Fill missing distance to lawn and arena edges
for j = 1:numTracks
    ind = data.wormNum == info.wormNum{plateInd}(j);
    data.distanceLawnEdge(ind) = fillmissing(data.distanceLawnEdge(ind),'linear',...
        'SamplePoints',time,'EndValues','nearest');
    data.distanceArenaEdge(ind) = fillmissing(data.distanceArenaEdge(ind),'linear',...
        'SamplePoints',time,'EndValues','nearest');
end

% Get identity and OD600 of closest lawn
mask = flip(info.lawnClosest{plateInd})';
data.closestLawnID = mask(indices);
mask = flip(info.lawnClosestOD600{plateInd})';
data.closestOD600 = mask(indices);

% Get tracks within 0.5 mm of arena boundary
data.nearArena = data.distanceArenaEdge./scale <= 0.5;

% Get tracks outside of boundary (untracked)
data.outOfBounds = false(size(data,1),1);
indCheck = find(data.nearArena & data.noTrack == 1);
for k = 1:length(indCheck)
    if data.nearArena(max(indCheck(k)-1,1)) == 1 || ...
            data.nearArena(min(indCheck(k)+1,size(data,1))) == 1 || ...
            data.outOfBounds(max(indCheck(k)-1,1)) == 1
        data.outOfBounds(indCheck(k)) = 1;
    end
end

end