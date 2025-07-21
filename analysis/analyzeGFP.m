function [] = analyzeGFP()
%[] = ANALYZEGFP()
%
%   ANALYZEGFP gets all the meta data for a set of experiments, identifies
%   brightfield images of acetate templates, identifies background images
%   of empty plates, and extracts and analyzes the fluorescence intensity
%   profile for each patch.  
%
%   Written 2/21/2024 by Jess Haley in MATLAB R2023b.
%
%   See also GETPLATEINFO, GETMETADATACZI, GETACETATEARENA, 
%   GETLAWNSTHRESHOLD, GETFLUORESCENCEBACKGROUND, ANALYZELAWNPROFILES.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Prepare Workspace and Get Experiments to Analyze

clear
close all

% Add paths
addpath(genpath('Z:\jhaley\code\bfmatlab'))
addpath(genpath('Z:\jhaley\foragingPaper'))

% Set info file
[infoFile,path] = uigetfile({'*.xls*','*.xlsx'},'Select experimentInfo file');

% Set directory to save all files in
saveDir = [path,'images',filesep];

% Get experiment info from notebook file
info = getPlateInfo([path,infoFile]);

% Get experiment ids
expToAnalyze = unique([info.expNum]);

% Query user to overwrite already saved data (if applicable)
overwrite = questdlg('Would you like to overwrite the saved data and metadata?',...
    'Overwrite','Yes','No','No');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Get MetaData and load images for the file(s)

metaData = [];
data = [];
fprintf('Loading metadata and images: \n')
for i = 1:length(expToAnalyze)
    tic

    % Get plate numbers for this experiment
    plateNums = info.plateNum([info.expNum] == expToAnalyze(i));

    % Get experiment date for this experiment
    expDate = datestr(mean(info.timeRoomTemp(plateNums),'omitnan'),'yy-mm-dd');
    
    % Build file paths
    expFolder = [saveDir,expDate,filesep];
    metaDataFile = [expFolder,'metadata.csv'];
    dataFile = [expFolder,'data.mat'];

    % Get metadata for each file
    if exist(metaDataFile,'file') && strcmp(overwrite,'No')
        thisMetaData = readtable(metaDataFile);
        thisMetaData.savedTime.TimeZone = 'local';
        thisMetaData.acquisitionTime.TimeZone = 'local';
        load(dataFile,'images');
    else
        [thisMetaData,images] = getMetaDataCZI(expFolder);
    end

    metaData = [metaData;thisMetaData];
    data = [data;images];
    
    fprintf('\t[Experiment %.2d: %.2d Images] %.1f s \n',...
        expToAnalyze(i),size(thisMetaData,1),toc)
end
clear images thisMetaData

% Add image number
metaData.imageNum(:) = 1:size(metaData,1);

% Add plate number (corresponds to folder name)
for i = 1:size(metaData,1)
    [~,subfolder,~] = fileparts(metaData.pathName(i));
    metaData.plateNum(i) = str2double(subfolder);
end

% Add info to metadata
metaData = join(metaData,info,'Keys','plateNum');

% Calculate growth times
metaData.growthTimeAfterSeed(:) = metaData.timeSeedColdRoom - ...
    metaData.timeSeed;
metaData.growthTimeColdRoom(:) = max(metaData.timeRoomTemp - ....
    metaData.timeSeedColdRoom,0);
metaData.timeRoomTemp.TimeZone = 'local';
metaData.growthTimeRoomTemp(:) = max(metaData.acquisitionTime - ...
    metaData.timeRoomTemp,0);
    
% Reorder metadata table
metaData = movevars(metaData,{'imageNum','plateNum','expNum'},...
    'Before','fileName');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Identify acetate templates from brightfield images

% Set template size
templateSize = [16,25];

% Initialize variables
template = struct('imageNum',[],'plateNum',[],'expNum',[],'mask',{[]},...
    'scale',[],'center',[],'offset',[],'angle',[]);
expScale = nan(size(expToAnalyze));

fprintf('Analyzing acetate templates: \n')
for i = 1:length(expToAnalyze)
    tic

    % Get experiment date for this experiment
    expDate = datestr(mean(metaData.savedTime(metaData.expNum == ...
        expToAnalyze(i)),'omitnan'),'yy-mm-dd');

    % Check if arena folder exists to save template
    arenaFolder = fullfile(saveDir,expDate,'data','arena');
    if ~exist(arenaFolder,'dir')
        mkdir(arenaFolder)
    end

    % Find brightfield images with acetate templates
    indTemplate = find(metaData.expNum == expToAnalyze(i) & ...
        ~metaData.fluorescence & strcmp(metaData.template,'rectangle'));

    for j = 1:length(indTemplate)
        ind = indTemplate(j);

        % Build file path to save analyzed arena
        arenaFile = fullfile(arenaFolder,...
            [metaData.fileName{ind}(1:end-3),'mat']);
        
        % Get template information
        if exist(arenaFile,'file') && strcmp(overwrite,'No')
            load(arenaFile,'arena');
        else
            arena = getAcetateArena(data{ind},templateSize);
            save(arenaFile,'arena');
        end
        
        % Add template info to template metadata
        template.imageNum(end+1) = metaData.imageNum(ind);
        template.plateNum(end+1) = metaData.plateNum(ind);
        template.expNum(end+1) = metaData.expNum(ind);
        template.mask{end+1} = arena.mask;
        template.scale(end+1) = arena.scale;
        template.center{end+1} = arena.center;
        template.offset{end+1} = arena.offset;
        template.angle(end+1) = arena.angle;
    end

    % Get average scale for each experiment
    expScale(i) = mean(template.scale(template.expNum == expToAnalyze(i)),'omitnan');
    metaData.scale(metaData.expNum == expToAnalyze(i)) = expScale(i);

    fprintf('\t[Experiment %.2d: %.2d Images] %.1f s \n',...
            expToAnalyze(i),length(indTemplate),toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Identify corresponding background and brightfield for each image

% Get matrix of durations between each image acquired
tAcquisition = seconds(duration(metaData.acquisitionTime - ...
    metaData.acquisitionTime(1)));
tDistance = squareform(pdist(tAcquisition));

% Get matrix of exposure times and plate numbers between each image acquired
exposureMatrix = squareform(pdist(metaData.exposureTime)) == 0;
plateMatrix = squareform(pdist(metaData.plateNum)) == 0;

% Find corresponding "background" images (i.e. image acquired with the 
% closest acquisition time, same exposure, and OD600 = 0)
indBackground = metaData.OD600 == 0 & strcmp(metaData.template,'none');
possibleBackground = tDistance;
possibleBackground(~exposureMatrix | ...
    repmat(~indBackground',size(indBackground))) = NaN;
[~,metaData.backgroundImageNum] = min(possibleBackground,[],2);

% Find corresponding brightfield image for each fluorescent image (i.e. 
% brightfield image of the same plate acquired with the closest acquisiton
% time; images without a template may not have a corresponding brightfield)
possibleBrightfield = tDistance;
possibleBrightfield(~plateMatrix | ...
    repmat(metaData.fluorescence',size(metaData.fluorescence)) | ...
    repmat(strcmp(metaData.template,'none')',size(metaData.fluorescence))) = NaN;
[~,metaData.brightfieldImageNum] = min(possibleBrightfield,[],2,'omitnan');
metaData.brightfieldImageNum(sum(~isnan(possibleBrightfield),2) == 0) = NaN;

% Save metadata to file
writetable(metaData,fullfile(saveDir,'metadata.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Get smoothed background

background = {};
backgroundInfo = table();
fprintf('Analyzing background images: \n')
for i = 1:length(expToAnalyze)
    tic

    % Get experiment date for this experiment
    expDate = datestr(mean(metaData.savedTime(metaData.expNum == ...
        expToAnalyze(i)),'omitnan'),'yy-mm-dd');

    % Check if background folder exists to save data
    backgroundFolder = fullfile(saveDir,expDate,'data','background');
    if ~exist(backgroundFolder,'dir')
        mkdir(backgroundFolder)
    end

    % Get background images for this experiment
    indBackground = unique(metaData.backgroundImageNum(...
        metaData.expNum == expToAnalyze(i) & metaData.fluorescence));

    for j = 1:length(indBackground)

        imageNum = indBackground(j);

        % Build file path to save analyzed background image
        backgroundFile = fullfile(backgroundFolder,...
            [metaData.fileName{imageNum}(1:end-3),'mat']);
        
        % Get template information
        if exist(backgroundFile,'file') && strcmp(overwrite,'No')
            load(backgroundFile,'imageFit','parameters');
        else
            [imageFit,parameters] = ...
                getFluorescenceBackground(data{imageNum},'smooth');
            save(backgroundFile,'imageFit','parameters');
        end

        background{end+1} = imageFit;
        backgroundInfo = [backgroundInfo;[table(imageNum),parameters]];
    end
    fprintf('\t[Experiment %.2d: %.2d Images] %.1f s \n',...
        expToAnalyze(i),length(indBackground),toc)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Extract and analyze edge profiles from flourescent images

v = VideoWriter([path,'fits'],'MPEG-4');
open(v);

lawnProfiles = struct('imageNum',[],'plateNum',[],'expNum',[],...
    'image',{[]},'imageNormalized',{[]},'labeled',{[]},...
    'distances',{[]},'pixelValues',{[]},'pixelValuesNormalized',{[]});
lawnAnalysis = table();

fprintf('Analyzing fluorescent images: \n')
for i = 1:length(expToAnalyze)
    tic

    % Get experiment date for this experiment
    expDate = datestr(mean(metaData.savedTime(metaData.expNum == ...
        expToAnalyze(i)),'omitnan'),'yy-mm-dd');

    % Check if arena folder exists to save metadata
    lawnFolder = fullfile(saveDir,expDate,'data','lawn');
    if ~exist(lawnFolder,'dir')
        mkdir(lawnFolder)
    end

    % Get indices of images to analyze
    indFluorescence = find(metaData.expNum == expToAnalyze(i) & ...
        metaData.fluorescence & metaData.OD600 > 0.1);

    for j = 1:length(indFluorescence)
        % Get image id
        imageNum = indFluorescence(j);

        % Build file path to save analyzed lawn
        lawnFile = fullfile(lawnFolder,...
            [metaData.fileName{imageNum}(1:end-3),'mat']);

        if exist(lawnFile,'file') && strcmp(overwrite,'No')
            load(lawnFile,'lawn');
        else
            % Get corresponding background image
            backgroundImage = background{backgroundInfo.imageNum == ...
                metaData.backgroundImageNum(imageNum)};

            % Get template (if applicable)
            indTemplate = find(template.imageNum == metaData.brightfieldImageNum(imageNum));
            if ~isempty(indTemplate)
                % Erode template mask by 0.1 mm
                templateMask = bwdist(~template.mask{indTemplate}) >= ...
                    0.1*metaData.scale(imageNum);
            else
                templateMask = true(size(backgroundImage));
            end

            % Get lawn radius and count
            lawnRadius = metaData.scale(imageNum)*...
                metaData.lawnDiameter(imageNum)/2;

            % Get lawn profiles
            lawn = analyzeLawnProfiles(data{imageNum},backgroundImage,...
                templateMask,lawnRadius,metaData.scale(imageNum));
        
            % Save lawn data
            save(lawnFile,'lawn');
        end

        % Add data to lawnProfiles
        lawnProfiles.imageNum(end+1) = metaData.imageNum(imageNum);
        lawnProfiles.plateNum(end+1) = metaData.plateNum(imageNum);
        lawnProfiles.expNum(end+1) = metaData.expNum(imageNum);
        lawnProfiles.image{end+1} = lawn.image;
        lawnProfiles.imageNormalized{end+1} = lawn.imageNormalized;
        lawnProfiles.labeled{end+1} = lawn.labeled;
        lawnProfiles.distances{end+1} = lawn.distances;
        lawnProfiles.pixelValuesNormalized{end+1} = lawn.pixelValuesNormalized;

        % Add analysis to lawnAnalysis
        if ~isempty(lawn.analysis)
            lawnAnalysis = [lawnAnalysis;...
                [table(imageNum*ones(size(lawn.analysis,1),1),...
                'VariableNames',{'imageNum'}),lawn.analysis]];
            figure(1); hold off;
            plot((lawn.distances - lawn.analysis.xOuterEdge)',...
                (lawn.pixelValuesNormalized - lawn.analysis.yOuterEdge)');
            hold on; scatter([lawn.analysis.xOuterEdge,lawn.analysis.xPeak]-lawn.analysis.xOuterEdge,...
                [lawn.analysis.yOuterEdge,lawn.analysis.yPeak] - lawn.analysis.yOuterEdge,'r')
            title(['imageNum = ',num2str(imageNum)])
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
    end

    fprintf('\t[Experiment %.2d: %.2d Images] %.1f s \n',...
        expToAnalyze(i),length(indFluorescence),toc)
end
writeVideo(v,frame);
close(v)

% Add metaData to lawnAnalysis
lawnAnalysis = join(lawnAnalysis,metaData,'key',{'imageNum'});
lawnAnalysis.growthTimeTotal = lawnAnalysis.growthTimeAfterSeed + ...
    lawnAnalysis.growthTimeRoomTemp;
lawnAnalysis.growthTimeCondition = 12*round(hours(lawnAnalysis.growthTimeTotal/12));
for i = 1:length(expToAnalyze)
    timeRoomTemp = max(info.timeRoomTemp(info.expNum == expToAnalyze(i)));
    timeRoomTemp.TimeZone = 'local';
    indExp = lawnAnalysis.expNum == expToAnalyze(i);
    lawnAnalysis.timeExperiment(indExp) = ...
        lawnAnalysis.acquisitionTime(indExp) - timeRoomTemp;
end

% Save workspace
save([path,'analyzeGFP_',char(datetime('today'),'yy-MM-dd'),'.mat'],'-v7.3')

end