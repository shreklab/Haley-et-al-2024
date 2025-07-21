% permutePatchLocation.m
%
%   permutePatchLocation generates a null model for foraging behavior by
%   randomizing the spatial locations of food patches. For each experimental
%   plate, the script creates thousands of new layouts by individually
%   rotating and translating the original patches to new, non-overlapping
%   locations within the arena. It then replays the worms' actual recorded
%   trajectories on these shuffled layouts to calculate a null distribution
%   of patch occupancy time, which can be used to assess whether the
%   observed behavior is significantly different from random chance.
%
%   CONFIGURATION:
%       - expName [char]: The name of the experiment to analyze, which
%           determines which data files are loaded (e.g.,
%           'foragingConcentration').
%
%   WORKFLOW:
%       1.  Load Data: Loads worm tracking data, experiment metadata
%           (including original patch masks), and encounter thresholds.
%       2.  Iterate Permutations: Loops through a set number of repetitions
%           (nReps).
%       3.  Shuffle Patches Per Plate: Within each repetition, it iterates
%           through each experimental plate. For each patch on that plate,
%           it randomly rotates and translates it to a new valid position
%           within the arena, ensuring no overlap with other shuffled patches.
%       4.  Calculate Null Occupancy: For each worm on that plate, it maps
%           its true trajectory onto the newly created shuffled patch layout
%           and interpolates the "on-patch" status over a one-hour period.
%       5.  Save Results: After all repetitions are complete, it saves the
%           resulting 'onPatchShuffled' logical matrix, which contains the
%           null occupancy data for every worm and every permutation.
%
%   REQUIRED DATA FILES:
%       - '<expName>/experimentInfo.mat': Metadata for the specified
%           experiment, including arena and lawn masks.
%       - '<expName>/midpoint.mat': Worm tracking data.
%       - 'foragingMini/encounterThresholds.csv': Contains thresholds for
%           defining patch boundaries.
%
%   OUTPUT FILES:
%       - '<expName>/permutePatchLocation.mat': A 3D matrix
%         (worm x time x repetition) containing the logical "on-patch"
%         status for each worm's real track on the shuffled patch layouts.
%
%   Written 7/11/2025 by Jess Haley in MATLAB R2024a.
%
%   See also ANALYZEFORAGING, GETFORAGINGINFO, IMTRANSLATE, IMROTATE, BWDIST.


% Load Data
expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
thresh = readtable([path,'foragingMini\encounterThresholds.csv']);

plateNums = unique(data.plateNum);
numPlates = length(plateNums);
wormNums = unique(data.wormNum);
numWorms = length(wormNums);

numSec = 3600;
oneHour = ((1:3600)-1)./60;
nReps = 1000;
% onPatchShuffled = false(numWorms,numSec,nReps);
noTrack = any([isnan(data.xPosition),isnan(data.yPosition)],2);
for n = 399:nReps
    for i = 1:numPlates
        % Get info for plate
        indInfo = info.plateNum == plateNums(i);
        lawnRadii = info.lawnRadii{indInfo};
        pixels = info.pixels{indInfo};
        arenaMask = info.arenaMask{indInfo};

        % Get x,y position of arena mask
        [xArena,yArena] = find(arenaMask);
        translateSpace = [xArena - pixels(1)/2,yArena - pixels(2)/2];

        % Shuffle patch locations
        lawnMaskShuffled = false(pixels);
        for j = 1:length(lawnRadii)
            % Get mask of current patch
            patchMask = info.lawnMask{indInfo}.*(info.lawnClosest{indInfo} == j);

            % Center patch, then randomly rotate and translate
            patchMaskCentered = imtranslate(patchMask,size(patchMask)./2 - info.lawnCenters{indInfo}(j,:),...
                'nearest','FillValues',0);
            patchMaskRotated = imrotate(patchMaskCentered,randi(360),'nearest','crop');
            patchMaskShuffled = imtranslate(patchMaskRotated,translateSpace(randi(length(translateSpace)),:),...
                'nearest','FillValues',0);

            % Check if translation overlaps w/ other patches and arena mask
            while sum((bwdist(patchMaskShuffled) <= -thresh.distMidpointEnter.*info.scale(indInfo)) & ...
                    ((bwdist(lawnMaskShuffled) <= -thresh.distMidpointEnter.*info.scale(indInfo)) | ...
                    ~arenaMask),'all') > 0
                patchMaskShuffled = imtranslate(patchMaskRotated,...
                    translateSpace(randi(length(translateSpace)),:),'nearest','FillValues',0);
            end
            lawnMaskShuffled = lawnMaskShuffled | patchMaskShuffled;
        end
        thisMaskShuffled = double((bwdist(lawnMaskShuffled) <= ...
            -thresh.distMidpointEnter.*info.scale(indInfo)));
        % figure; imshowpair(info.lawnMask{indInfo} | ~arenaMask,...
        %     lawnMaskShuffled + (bwdist(lawnMaskShuffled) <= -thresh.distMidpointEnter.*info.scale(indInfo)))
        
        plateWorms = unique(data.wormNum(data.plateNum == plateNums(i)));
        for j = 1:length(plateWorms)

            % Get index
            indData = data.wormNum == plateWorms(j) & ~noTrack;
            indPixel = sub2ind(pixels,round(data.yPosition(indData)),...
                round(data.xPosition(indData)));
            onPatchShuffled(wormNums == plateWorms(j),:,n) = interp1(data.timeOffset(indData),...
                thisMaskShuffled(indPixel),oneHour.*60) >= 0.5;
        end
    end
end
save(fullfile(path,expName,'permutePatchLocation.mat'),'onPatchShuffled','-v7.3');