function [] = analyzeForaging()
%[] = ANALYZEFORAGING()
%
%   ANALYZEFORAGING uses the body segment location of an animal combined
%   with the location of bacterial patches/lawns to identify encounters and 
%   properties of those encounters.
%
%   INPUTS:
%       - infoFile: a .mat file that includes metadata for each
%           experiment. See GETFORAGINGINFO for more information on
%           expected variables.
%
%   FUNCTIONS:
%       ANALYZE BEHAVIOR
%       - analyzeWormLabTracks: takes body segment data exported from 
%           WormLab and computes numerous metrics related to the worm's 
%           position, velocity, and turning behavior as well as it's
%           location relative to the arena and lawns.
%       - analyzeEncounters: uses the location of animal(s) and lawn(s) to 
%           identify encounters and properties of those encounters.
%
%       CREATE PLOTS
%       - plotTracks: takes in info [struct] and data [struct] and plots 
%           the tracks of the specified wormNum [int] colored by the given
%           metric ('time', 'velocity', or 'path angle') [str] displayed on
%           top of an image of the patches and arena; saved as .pdf file to
%           the specified saveDir [str]
%       - createVideo: writes a video to file showing either the tracks of 
%           an animal or the original video downsampled with scale bar and
%           time stamp showing.
%
%   Written 1/11/2023 by Jess Haley in MATLAB R2023a.
%
%   See also GETFORAGINGINFO, ANALYZEWORMLABTRACKS, OFFSETTIME,
%   ANALYZEENCOUNTERS, DEFINEENCOUNTER, PLOTTRACKS, CREATEVIDEO.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Prepare Workspace and Get Experiments to Analyze

% Add paths
addpath(genpath('F:\'))
addpath(genpath('Z:\jhaley\foragingPaper\'))

% Set info file
[infoFile,path] = uigetfile('*.mat*');
load([path,infoFile],'info')
[~,expName] = fileparts(path(1:end - ...
        double(find(path == filesep,1,'last') == length(path))));

% Set directory to save all files in
saveDir = [path,'videos',filesep];

% Query user for experiment number(s) (if not specified as an input)
[expToAnalyze] = listdlg('PromptString',{'Select experiment(s) to analyze.'},...
    'ListString',num2str(unique([info.expNum])));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Load and Analyze X,Y Positions from Exported .xls Files from WormLab

overwrite = questdlg('Would you like to overwrite the Worm Lab data?',...
    'Overwrite','Yes','No','No');

bodyPart = 'midpoint';
if strcmp(expName,'foragingMini')
    bodyPart = questdlg('Would you like to analyze the head or midpoint?',...
        'Body Part','head','midpoint','tail','midpoint');
end

data = table();
fprintf('Analyzing WormLab Data: \n')
for i = 1:length(expToAnalyze)
    
    tic
    fprintf('\t[Experiment %.2d]: [Plate',expToAnalyze(i))
    
    % Get plate numbers for this experiment
    plateNums = info.plateNum([info.expNum] == expToAnalyze(i));
    videoNums = info.videoNum([info.expNum] == expToAnalyze(i));
    numPlates = length(plateNums); % # of plates in this experiment
    expDate = datestr(mean(info.timeRecord(info.expNum == expToAnalyze(i)),...
        'omitnan'),'yy-mm-dd');

    % Check if track folder exists to save data
    wormlabFolder = fullfile(saveDir,expDate,'wormlab',filesep);
    if strcmp(expName,'foragingMini')
        wormlabFolder = fullfile(wormlabFolder,'centerpoints',filesep);
    end
    trackFolder = fullfile(saveDir,expDate,'data',bodyPart);
    if ~exist(trackFolder,'dir')
        mkdir(trackFolder)
    end

    for j = 1:numPlates
        plateNum = plateNums(j);
        videoNum = videoNums(j);
        ind = info.plateNum == plateNum & info.videoNum == videoNum;

        if isempty(info.videoFileName{ind})
            continue
        else
            % Build file path to save analyzed tracks
            trackFile = fullfile(trackFolder,[info.videoFileName{ind}(1:end-3),'csv']);

            % Analyze tracks
            if exist(trackFile,'file') && strcmp(overwrite,'No')
                newData = readtable(trackFile);
            elseif exist([wormlabFolder,info.wormLabFileName{ind}],'file') == 0
                newData = table();
            else
                newData = analyzeWormLabTracks(info,ind,wormlabFolder,bodyPart);
                writetable(newData,trackFile);
            end

            % Append data to table
            data = [data;newData];
            fprintf(' %.3d',plateNum)
        end
    end

    fprintf('] %.1f s \n',toc)
end

% Merge worm tracks from multi-hour/multi-video recordings
data = offsetTime(data,info);

save([path,bodyPart,'.mat'],'data','-v7.3');
fprintf('Data saved: %s \n',[path,bodyPart,'.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Analyze Encounters

overwrite = questdlg('Would you like to overwrite the analyses?',...
    'Overwrite','Yes','No','No');

encounter = table(); % rows = # lawn visits
distance = struct('info',[],'probReside',[],'velocitySmooth',[]); % columns = # distance bins
trajectory = struct('enterInfo',[],'enterVelocity',[],'enterAngle',[],'enterDistLawn',[],'enterDistEnter',[],...
    'exitInfo',[],'exitVelocity',[]); % timepoints immediately before/after lawn encounter

fprintf('Analyzing Behavior: \n')
for i = 1:length(expToAnalyze)
    
    tic
    fprintf('\t[Experiment %.2d]: [Plate',expToAnalyze(i))
    
    % Get plate numbers for this experiment
    plateNums = unique(info.plateNum([info.expNum] == expToAnalyze(i)));
    numPlates = length(plateNums); % # of plates in this experiment
    expDate = datestr(mean(info.timeRecord(info.expNum == expToAnalyze(i)),...
        'omitnan'),'yy-mm-dd');

    % Check if analysis folder exists to save metadata
    analysisFolder = fullfile(saveDir,expDate,'analysis',bodyPart);
    if ~exist(analysisFolder,'dir')
        mkdir(analysisFolder)
    end

    % Check if encounter folder exists to save metadata
    encounterFolder = fullfile(analysisFolder,'encounter');
    if ~exist(encounterFolder,'dir')
        mkdir(encounterFolder)
    end

    % Check if distance folder exists to save metadata
    distanceFolder = fullfile(analysisFolder,'distance');
    if ~exist(distanceFolder,'dir')
        mkdir(distanceFolder)
    end

    % Check if border folder exists to save metadata
    trajectoryFolder = fullfile(analysisFolder,'trajectory');
    if ~exist(trajectoryFolder,'dir')
        mkdir(trajectoryFolder)
    end

    for j = 1:numPlates
        plateNum = plateNums(j);
        if sum(data.plateNum == plateNum) == 0
            continue
        else
        
        ind = info.plateNum == plateNum & info.videoNum == 1;
            
        % Build file path to save analyzed data
        encounterFile = fullfile(encounterFolder,...
            [info.videoFileName{ind}(1:end-3),'csv']);
        distanceFile = fullfile(distanceFolder,...
            [info.videoFileName{ind}(1:end-3),'mat']);
        trajectoryFile = fullfile(trajectoryFolder,...
            [info.videoFileName{ind}(1:end-3),'mat']);
        
        % Analyze tracks
        if exist(encounterFile,'file') && strcmp(overwrite,'No')
            encounterNew = readtable(encounterFile);
            distanceNew = load(distanceFile);
            trajectoryNew = load(trajectoryFile);
        else
            [encounterNew,distanceNew,trajectoryNew] = analyzeEncounters(info,data,plateNum);
            writetable(encounterNew,encounterFile);
            save(distanceFile,'-struct','distanceNew');
            save(trajectoryFile,'-struct','trajectoryNew');
        end

        % Append data to tables/structs
        encounter = [encounter; encounterNew];
        vars = fieldnames(distance);
        for j = 1:length(vars)
            distance.(vars{j}) = [distance.(vars{j});distanceNew.(vars{j})];
        end
        vars = fieldnames(trajectory);
        for k = 1:length(vars)
            trajectory.(vars{k}) = [trajectory.(vars{k});trajectoryNew.(vars{k})];
        end
        fprintf(' %.3d',plateNum)
        end
    end

    fprintf('] %.1f s \n',toc)
end

% Add more data to encounters
encounter.expName(:) = {expName};
encounter = join(encounter,info(:,{'plateNum','videoNum',...
        'exclude','strainName','strainID','peptone','OD600Label','lawnVolume',...
        'growthLawnGrowth','wormGrowth'}),'Keys',{'plateNum','videoNum'});
encounter.growthLawnGrowth = encounter.growthLawnGrowth + encounter.wormGrowth;
encounter.growthLawnGrowth(isnan(encounter.growthLawnGrowth)) = ...
    mean(encounter.growthLawnGrowth,'omitnan');

% Save encounters
save([path,'encounter.mat'],'encounter','-v7.3');
save([path,'distance.mat'],'distance','-v7.3');
save([path,'trajectory.mat'],'trajectory','-v7.3');
fprintf('Encounters saved: %s \n',[path,'encounter.mat'])

clear analysisNew newData expDate i j k numPlates plateNum plateNums vars;
save([path,'analyzeForagingWorkspace_',char(datetime('today'),'yy-MM-dd'),'.mat'],'-v7.3')

%% Plot paths

% Plot paths colored by time, velocity, and path angle
wormNums = unique(data.wormNum);
tic
for i = 1:length(wormNums)
    plotTracks(info,data,wormNums(i),'timeOffset',saveDir,bodyPart,overwrite);
    plotTracks(info,data,wormNums(i),'velocitySmooth',saveDir,bodyPart,overwrite);
    plotTracks(info,data,wormNums(i),'pathAngle',saveDir,bodyPart,overwrite);
    fprintf(' %.4d',wormNums(i))
    if mod(i,10)==0
        fprintf(' - %.1f s \n',toc)
        tic
    end
end

%% Create Video of worm tracks

createVideo(data,info);

end