function [encounter,distance,trajectory] = analyzeEncounters(info,data,plateNum)
%[encounter,distance,trajectory] = ANALYZEENCOUNTERS(info,data,plateNum)
%
%   ANALYZEENCOUNTERS uses the location of animal(s) and lawn(s) on a
%   given plate to identify encounters and properties of those encounters.
%
%   INPUTS:
%       - info [struct]: structure containing the metadata for all worm(s);
%           created with GETFORAGINGINFO
%       - data [struct]: structure containing the location of worm(s) at 
%           every timepoint; created with ANALYZEFORAGING
%       - plateNum [double]: identifying plate # to analyze
%
%   OUTPUTS:
%       - encounter [table]: a table containing the following variables
%           - expNum [uint16]: unique number identifier of each 
%               experimental date; plates on this date share the same
%               growth conditions (e.g. 1)
%           - plateNum [uint16]: unique number identifier of each 
%               experimental plate; worms on these plates share the same
%               growth and experimental conditions (e.g. 1)
%           - videoNum [uint16]: unique number identifier for each video
%               recorded for each plate (e.g. 1); only applicable for
%               experiments where multiple videos were recorded
%           - wormNum [uint16]: unique number identifier for each worm
%               (e.g. 1)
%           - id [double]: unique number identifier for each encounter for
%               each worm (e.g. 1); 0 indicates off-patch metrics before
%               the first encounter
%           - enter [double]: first frame of the encounter
%           - exit [double]: last frame of the encounter
%           - timeEnter [double]: time that the encounter began in seconds
%           - timeExit [double]: time that the encounter ended in seconds
%           - duration [double]: duration of the encounter (i.e. timeExit -
%               timeEnter) in seconds
%           - lawnID [double]: unique nuymber identifier of the bacterial
%               lawn visited during the encounter; matches
%               info.lawnClosest (e.g. 1)
%           - lawnOD600 [double]: OD600 pipetted on the plate; matches
%               info.lawnClosestOD600 (e.g. 1, 5, or 10)
%           - growthCondition [double]: intended hours of growth (i.e. 0, 
%               12, or 48)
%           - lawnGrowth [double]: duration of time lawns were allowed to
%               grow for at room temperature after drying in seconds;
%               matches the time for the start of the experiment
%           - velocityOn [double]: median velocity of animal on patch
%               during encounter
%           - velocityOnMin [double]: minimum velocity of animal on patch
%               during encounter
%           - distanceOnMax [double]: maximum distance animal's midpoint
%               travelled off of patch during the encounter; negative
%               values represent distance from patch edge in mm
%           - velocityBeforeEnter [double]: peak velocity within 10 seconds
%               before the start of the encounter
%           - velocityAfterEnter [double]: minimum velocity within 10
%               seconds after the start of the encounter
%           - timeSlowDown [double]: duration of time between the peak
%               velocity before and minimum velocity after the start of the
%               encounter in seconds
%           - decelerate [double]: slope of line fit to the velocity from
%               -1.5 seconds before to +6.5 seconds after the start of the
%               encounter
%           - velocityOff [double]:  median velocity of animal in between
%               encounters
%           - distanceOffMax [double]: maximum distance animal's midpoint
%               travelled from the patch edge in mm
%           - censorEnter [logical]: true for encounters that started on or 
%               before the first frame
%           - censorExit [logical]: true for encounters that ended on or 
%               after the last frame
%           - expName [str]: unique identifier for the name of the
%               experiment
%           - exclude [logical]: true if experiments are to be excluded
%               from future analyses (e.g. experiments lacking a contrast
%               video where lawns cannot be reliably detected)
%           - strainName [str]: unique identifier for the name of the
%               strain (e.g. well-fed, food-deprived, or osm-6)
%           - strainID [str]: unique identifier for the CGC strain id (e.g.
%               N2, or PR811)
%           - peptone [str]: describes whether the condition plate was made
%               with or without peptone (i.e. 'with' or 'without')
%           - OD600Label [str]: OD600 solution pipetted on each plate (e.g.
%               '1.00' or '0.50'); for experiments with multiple patch 
%               concentrations, solutions are listed (e.g. '1.00 5.00
%               10.00')
%           - lawnVolume [double]: pipetted volume of each lawn in uL 
%               (e.g. 0.5)
%           - growthlawnGrowth [duration]: duration of time growth plate
%               lawns were allowed to grow for at room temperature after
%               drying (i.e. HH:MM:SS)
%           - wormGrowth [duration]: duration of time between L4s being
%               picked and the experiment start time (i.e. HH:MM:SS)
%       - distance [struct]: a structure containing the following fields
%           - info [table]: a table containing identifying information for
%               each worm (i.e. expNum, plateNum, videoNum, wormNum) 
%               corresponding to the data in probReside and velocitySmooth
%           - probReside [Wx371 double]: the fraction of time that each of 
%               W worms resides within a bin of specified distance from the
%               patch edge; bins are span the array [-14.025:0.05:4.525] mm
%           - velocitySmooth [Wx371 double]: the mean velocity that each of 
%               W worms maintained within a bin of specified distance from 
%               the patch edge; bins are span the array 
%               [-14.025:0.05:4.525] mm
%       - trajectory [struct]: a structure containing the following fields
%           - enterInfo [table]: a table containing identifying information
%               for each patch enter event (start of encounter) (i.e.
%               expNum, plateNum, videoNum, wormNum, id) corresponding to
%               the data in enterVelocity, enterAngle, enterDistLawn, and
%               enterDistEnter
%           - enterVelocity [Ex1201 double]: interpolated velocity for a
%               time window (i.e. [-60:0.25:240] seconds) surrounding the
%               start of each of E patch entry events 
%           - enterAngle [Ex1201 double]: interpolated path angle for a
%               time window (i.e. [-60:0.25:240] seconds) surrounding the
%               start of each of E patch entry events
%           - enterDistLawn [Ex1201 double]: interpolated distance between 
%               the animal's midpoint and the closes point along the patch
%               edge for a time window (i.e. [-60:0.25:240] seconds)
%               surrounding the start of each of E patch entry events
%           - enterDistEnter [Ex1201 double]: interpolated distance between 
%               the animal's midpoint and the location of patch entry for a
%               time window (i.e. [-60:0.25:240] seconds) surrounding the
%               start of each of E patch entry events
%           - exitInfo [table]: a table containing identifying information
%               for each patch exit event (end of encounter) (i.e.
%               expNum, plateNum, videoNum, wormNum, id) corresponding to
%               the data in exitVelocity
%           - exitVelocity [Ex1201 double]: interpolated velocity for a
%               time window (i.e. [-60:0.25:240] seconds) surrounding the
%               end of each of E patch exit events 
%
%   Written 3/26/2024 by Jess Haley in MATLAB R2024a.
%
%   See also ANALYZEFORAGING, GETFORAGINGINFO, DEFINEENCOUNTER, CSAPS,
%   GRADIENT, MOVMIN, INTERP1, HISTCOUNTS.

warning('off','SPLINES:CHCKXYWP:NaNs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Initialize Analysis Table with Identifying Information and Metadata

wormNums = info.wormNum{[info.plateNum] == plateNum};
numWorms = length(wormNums); % # of worms on this plate

% Get thresholds
thresh = readtable('Z:\jhaley\foragingPaper\foragingMini\encounterThresholds.csv');

% Initialize table
encounter = table(); % rows = # lawn visits
distance = struct('info',[],'probReside',[],'velocitySmooth',[]); % columns = # distance bins
trajectory = struct('enterInfo',[],'enterVelocity',[],'enterAngle',[],'enterDistLawn',[],'enterDistEnter',[],...
    'exitInfo',[],'exitVelocity',[]); % timepoints immediately before/after lawn encounter

% Get scale and frame rate
indInfo = find(info.plateNum == plateNum);
scale = sum(info.scale(indInfo).*info.numFrames(indInfo))./...
    sum(info.numFrames(indInfo)); % weighted mean if more than 1 video
frameRate = sum(info.frameRate(indInfo).*info.numFrames(indInfo))./...
    sum(info.numFrames(indInfo)); % weighted mean if more than 1 video

for i = 1:numWorms

    % Get indices for this worm in the data table
    wormNum = wormNums(i);
    ind = find(data.wormNum == wormNum);

    if isempty(ind)
        continue
    end

    % Initialize tables
    encounterWorm = table();
    trajectoryWorm = struct();
    distanceWorm = struct();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Get Lawn Encounters (i.e. Enter/Exit Events)

    % Get lawn IDs from data
    videoNum = data.videoNum(ind);
    timeOffset = data.timeOffset(ind);

    % Get data for worm
    distanceLawnEdge = data.distanceLawnEdge(ind)./scale;
    closestID = data.closestLawnID(ind);
    closestOD600 = data.closestOD600(ind);
    velocitySmooth = data.velocitySmooth(ind);
    pathAngle = data.pathAngle(ind);

    % Smooth the velocity using cubic spline interpolation and
    % compute gradient
    velocitySpline = csaps(timeOffset,velocitySmooth,0.5,timeOffset);
    accelerationSpline = gradient(velocitySpline,timeOffset);
    
    if strcmp(data.bodyPart(ind(1)),'head')
        % Get enter and exit events where head is on lawn
        bodyPartEnter = find(diff([0;distanceLawnEdge >= 0])==1);
        bodyPartExit = find(diff([distanceLawnEdge >= 0;0])==-1);

    else
        % Get enter and exit events where head should be able to touch
        % lawn (~0.45 mm from patch edge)
        bodyPartEnter = find(diff([0;distanceLawnEdge >= thresh.distMidpointEnter])==1);
        bodyPartExit = find(diff([distanceLawnEdge >= thresh.distMidpointEnter;0])==-1);

        % Remove events where midpoint never gets close enough to lawn
        % (~0.3 mm from patch edge)
        if ~isempty(bodyPartEnter)
            fullEnter = arrayfun(@(ent1,ent2) ...
                sum(distanceLawnEdge(ent1:ent2) >= thresh.distMidpointMin),...
                bodyPartEnter,[bodyPartEnter(2:end);length(ind)]);
            bodyPartEnter = bodyPartEnter(fullEnter > 0);
        end
    end

    % Iteratively add enter and exit events
    if ~isempty(bodyPartEnter) & ~isempty(bodyPartExit)
        if bodyPartEnter(1) == 1
            wormEnter = 1;
            wormExit = bodyPartExit(1);
        else
            wormEnter = 1;
            wormExit = 1;
        end
        while wormExit(end) < bodyPartExit(end) && wormExit(end) < bodyPartEnter(end)
            wormEnter(end+1,1) = bodyPartEnter(find(bodyPartEnter > wormExit(end),1,'first'));
            wormExit(end+1,1) = bodyPartExit(find(bodyPartExit >= wormEnter(end),1,'first'));
        end

        % Check distance variability off patch & remove exit events
        % matching on patch variability
        distanceVariability = arrayfun(@(ent,ex) ...
            std(distanceLawnEdge(ex:ent)),wormEnter(2:end),wormExit(1:end-1));
        removeExit = find(distanceVariability < thresh.distVarMax);
        lawnID = reshape(closestID([wormEnter(removeExit+1),wormExit(removeExit)]),[],2);
        removeExit(lawnID(:,1)~=lawnID(:,2)) = [];
        wormEnter(removeExit+1) = [];
        wormExit(removeExit) = [];
    else
        wormEnter = 1;
        wormExit = 1;
    end
    encounterWorm.enter = wormEnter;
    encounterWorm.exit = wormExit;

    % Calculate timing of encounter
    encounterWorm.timeEnter = timeOffset(encounterWorm.enter);
    encounterWorm.timeExit = timeOffset(encounterWorm.exit);
    encounterWorm.duration(:) = encounterWorm.timeExit(:) - ...
        encounterWorm.timeEnter(:);

    % Get info about encounters
    encounterWorm.lawnID = closestID(encounterWorm.enter);
    encounterWorm.lawnOD600 = closestOD600(encounterWorm.enter);
    encounterWorm.growthCondition(:) = info.growthCondition(indInfo(1));
    encounterWorm.lawnGrowth = minutes(info.lawnGrowth(indInfo(1))) + encounterWorm.timeEnter./60;

    % Remove the middle lawn (id = 10) if concentration is 0
    encounterWorm.lawnOD600(isnan(encounterWorm.lawnID)) = NaN;
    if max(info.lawnClosestOD600{indInfo(1)},[],'all') > 1e-8
        encounterWorm(encounterWorm.lawnOD600 == 0,:) = [];
    end

    % Get velocity on patch and minimum distance
    winSize = round(10*frameRate);
    velocityAfter = movmin(velocitySmooth,[0 winSize]);
    for k = 1:height(encounterWorm)
        encounterWorm.velocityOn(k) = ...
            median(velocitySmooth(encounterWorm.enter(k):encounterWorm.exit(k)),'omitnan');
        encounterWorm.velocityOnMin(k) = ...
            min(velocitySmooth(encounterWorm.enter(k):encounterWorm.exit(k)));
        encounterWorm.distanceOnMax(k) = ...
            max(distanceLawnEdge(encounterWorm.enter(k):encounterWorm.exit(k)));
        
        % Get max velocity at peak before enter
        peakBefore = find(accelerationSpline >= 0 & timeOffset <= encounterWorm.timeEnter(k) & ...
            timeOffset >= encounterWorm.timeEnter(k) - 10,1,'last');
        if isempty(peakBefore)
            peakBefore = encounterWorm.enter(k);
        end
        peakBefore = find(abs(timeOffset - timeOffset(peakBefore)) <= 1);
        [encounterWorm.velocityBeforeEnter(k),indBefore] = max(velocitySmooth(peakBefore));
        
        % Get min velocity after enter
        encounterWorm.velocityAfterEnter(k) = max([velocityAfter(encounterWorm.enter(k)),...
            encounterWorm.velocityOnMin(k)],[],2);
        indAfter = max([1,find(velocitySmooth == encounterWorm.velocityAfterEnter(k) & ...
            timeOffset >= encounterWorm.timeEnter(k) & timeOffset <= encounterWorm.timeExit(k),1,'first')]);

        % Duration of slowdown
        encounterWorm.timeSlowDown(k) = timeOffset(indAfter) - ...
            timeOffset(peakBefore(indBefore));

        % Compute deceleration at slowdown via linear fit from -1.5 to +6.5
        % seconds from patch entry
        indDecelerate = timeOffset >= encounterWorm.timeEnter(k) - 1.5 & ...
            timeOffset <= encounterWorm.timeEnter(k) + 6.5;
        dT = timeOffset(indDecelerate) - mean(timeOffset(indDecelerate),'omitnan');
        dV = velocitySmooth(indDecelerate) - mean(velocitySmooth(indDecelerate),'omitnan');
        encounterWorm.decelerate(k) = sum(dT.*dV,'omitnan')./sum(dT.^2,'omitnan');
    end

    % Get velocity off patch and maximum distance worm travels off lawn between events
    encounterWorm.velocityOff(:) = NaN;
    encounterWorm.distanceOffMax(:) = NaN;
    for k = 1:height(encounterWorm) - 1
        encounterWorm.velocityOff(k) = ...
            median(velocitySmooth(encounterWorm.exit(k):encounterWorm.enter(k+1)),'omitnan');
        encounterWorm.distanceOffMax(k) = ...
            abs(min(distanceLawnEdge(encounterWorm.exit(k):encounterWorm.enter(k+1))));
    end

    % Add experiment identifiers
    encounterWorm.expNum(:) = info.expNum(indInfo(1));
    encounterWorm.plateNum(:) = plateNum;
    encounterWorm.videoNum = videoNum(encounterWorm.enter);
    encounterWorm.wormNum(:) = wormNum;
    
    % Off data in the first row corresponds to the time before the first
    % encounter; on data is not relevant if exit == 1
    if encounterWorm.exit(1) == 1
        resetVars = {'enter','exit','timeEnter','timeExit','duration',...
            'lawnID','lawnOD600','growthCondition','lawnGrowth',...
            'velocityOn','velocityOnMin','distanceOnMax',...
            'velocityBeforeEnter','velocityAfterEnter','timeSlowDown','decelerate'};
        encounterWorm(1,resetVars) = num2cell(nan(size(resetVars)));
    end

    % Off data in the last row corresponds to the time after the last encounter
    % (if applicable)
    if encounterWorm.exit(end) == length(ind) % if the last event contains the last frame, remove off data
        resetVars = {'velocityOff','distanceOffMax'};
        encounterWorm(end,resetVars) = num2cell(nan(size(resetVars)));
    end

    % Label encounters that are censored
    encounterWorm.censorEnter = encounterWorm.enter == 1; % left-censored (timeEnter <= 0)
    encounterWorm.censorExit = encounterWorm.exit == length(ind); % right-censored (max(timeOffset) <= timeExit)

    % Reorder variables
    encounterWorm.id = cumsum(encounterWorm.enter > 0);
    encounterWorm = movevars(encounterWorm,{'expNum','plateNum',...
        'videoNum','wormNum','id'},'Before','enter');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. Add Time-Series Analyses about Behavior upon Enter/Exit Lawn

    % Define time window
    winTime = [60 240]; % -60 to +180 seconds from enter/exit event
    timeWindow = -winTime(1):0.25:winTime(2);

    % Get indices of time points in window surrounding enter/exit events
    includeEnter = encounterWorm.enter > 0; % do I want this to be 1?
    includeExit = encounterWorm.exit < length(ind);
    timeEnter = encounterWorm.timeEnter(includeEnter);
    timeExit = encounterWorm.timeExit(includeExit);

    % Get trajectories of entrances
    trajectoryWorm.enterInfo = encounterWorm(includeEnter,...
        {'expNum','plateNum','videoNum','wormNum','id'});
    enterVelocity = arrayfun(@(t) interp1(timeOffset,velocitySmooth,t + timeWindow),...
        timeEnter,'UniformOutput',false);
    trajectoryWorm.enterVelocity = vertcat(enterVelocity{:});
    enterAngle = arrayfun(@(t) interp1(timeOffset,pathAngle,t + timeWindow),...
        timeEnter,'UniformOutput',false);
    trajectoryWorm.enterAngle = vertcat(enterAngle{:});
    enterDistLawn = arrayfun(@(t) interp1(timeOffset,distanceLawnEdge,t + timeWindow),...
        timeEnter,'UniformOutput',false);
    trajectoryWorm.enterDistLawn = vertcat(enterDistLawn{:});
    trajectoryWorm.enterDistEnter = nan(length(timeEnter),length(timeWindow));
    for k = 1:length(timeEnter)
        xPosition = interp1(timeOffset,data.xPosition(ind),timeEnter(k) + timeWindow);
        yPosition = interp1(timeOffset,data.yPosition(ind),timeEnter(k) + timeWindow);
        trajectoryWorm.enterDistEnter(k,:) = sqrt((xPosition - xPosition(timeWindow == 0)).^2 + ...
            (yPosition - yPosition(timeWindow == 0)).^2)./scale;
    end

    % Get trajectories of exits
    trajectoryWorm.exitInfo = encounterWorm(includeExit,...
        {'expNum','plateNum','videoNum','wormNum','id'});
    exitVelocity = arrayfun(@(t) interp1(timeOffset,velocitySmooth,t + timeWindow),...
        timeExit,'UniformOutput',false);
    trajectoryWorm.exitVelocity = vertcat(exitVelocity{:});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. Add Distance Analyses about Average Behavior near/on Lawn

    % Define distance window
    binEdges = -14.025:0.05:4.525; % mm
    % distanceWorm.distLawnEdge = binEdges(1:end-1)+0.05;

    % Bin data by distance from lawn edge
    [binN,~,binID] = histcounts(distanceLawnEdge,binEdges);
    binExist = binN > 0;
    G = findgroups(binID);

    % Get bin metrics
    distanceWorm.info = encounterWorm(1,{'expNum','plateNum','videoNum','wormNum'});
    distanceWorm.probReside = binN./length(ind);
    distanceWorm.velocitySmooth = nan(size(distanceWorm.probReside));
    distanceWorm.velocitySmooth(binExist) = splitapply(@(x) mean(x,'omitnan'), velocitySmooth, G);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. Concatenate Analyses

    encounter = [encounter; encounterWorm];
    vars = fieldnames(distance);
    for j = 1:length(vars)
        distance.(vars{j}) = [distance.(vars{j});distanceWorm.(vars{j})];
    end
    vars = fieldnames(trajectory);
    for j = 1:length(vars)
        trajectory.(vars{j}) = [trajectory.(vars{j});trajectoryWorm.(vars{j})];
    end
end

end