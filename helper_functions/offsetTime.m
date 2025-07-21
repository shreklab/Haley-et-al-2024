function [data] = offsetTime(data,info)
%[data] = OFFSETTIME(data,info)
%
%   OFFSETTIME corrects the wormNum identifiers and time values for 
%   experiments where the recording spans multiple .avi files.
%
%   Written 5/8/2023 by Jess Haley in MATLAB R2023a.
%
%   See also ANALYZEFORAGING, GETEXPERIMENTINFO, GETFORAGINGNFO.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Get Time Offset

timeOffsets = zeros(size(data,1),1);
wormNums = unique(data.wormNum);
for i = 1:length(wormNums)
    
    % Get video numbers for this worm
    plateNum = unique(data.plateNum(data.wormNum == wormNums(i)));
    videoNums = unique(info.videoNum(info.plateNum == plateNum));

    for j = 1:length(videoNums)

        % Get time offsets (i.e. time between start of this recording and 
        % time of first recording containing this wormNum)
        ind = data.wormNum == wormNums(i) & data.videoNum == videoNums(j);
        timeOffsets(ind) = seconds(info.timeRecord(info.plateNum == plateNum & ...
            info.videoNum == videoNums(j)) - min(info.timeRecord(info.plateNum == plateNum)));
    end
end

% Add offsets to data
timeOffset = data.time + timeOffsets;
data = addvars(data,timeOffset,'After','time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Renumber wormNums to match across recordings

% Get ids for plates with time offsets (not the first plate)
plateNums = info.plateNum(info.videoNum > 1);
videoNums = info.videoNum(info.videoNum > 1);
for j = 1:length(plateNums)

    % Get worm and plate ids that need to be matched
    theseWormNums = unique(data.wormNum(data.plateNum == plateNums(j)));

    % Find index of first tracked position for each worm (some worms may
    % not have been tracked in the first frame(s) of the video)
    indFirstFrames = arrayfun(@(wormNum) find(data.wormNum == wormNum & ...
            data.plateNum == plateNums(j) & data.videoNum == videoNums(j) & ...
            ~isnan(data.xPosition),1,'first'),theseWormNums);

    % Find index of last tracked position for each worm in the last video
    try
    indLastFrames = arrayfun(@(wormNum) find(data.wormNum == wormNum & ...
            data.plateNum == plateNums(j) & data.videoNum == videoNums(j)-1 & ...
            ~isnan(data.xPosition),1,'last'),theseWormNums);
    catch ME
        if strcmp(ME.identifier,'MATLAB:arrayfun:NotAScalarOutput')
            fprintf('Time offsets not found for plate # %.3d \n',plateNums(j))
            continue
        end
    end

    % Get distances between each first and last frame worm
    distFrames = pdist2([data.xPosition(indFirstFrames),data.yPosition(indFirstFrames)],...
        [data.xPosition(indLastFrames),data.yPosition(indLastFrames)]);

    % Get all possible combinations of reordering the worm numbers
    combos = perms(1:length(theseWormNums));

    % Get distances of all possible combinations
    indCombos = sub2ind(size(distFrames),repmat(1:4,size(combos,1),1),combos);
    distCombos = distFrames(indCombos);

    % Best combo is the one with the smallest sum of distances
    [~,bestCombo] = min(sum(distCombos,2));

    % Reassign the wormNum ids
    newWormNum = theseWormNums(combos(bestCombo,:));
    theseWormInds = arrayfun(@(wormNum) find(data.wormNum == wormNum & ...
            data.videoNum == videoNums(j)),theseWormNums,'UniformOutput',false);
    for i = 1:length(theseWormNums)
        data.wormNum(theseWormInds{i}) = newWormNum(i);
    end
end

end