function [info] = getForagingInfo()
%[info] = GETFORAGINGINFO()
%
%   GETFORAGINGINFO compiles all of the metadata for a set of experiments
%   outlined in an infoFile (.xls) and saves that info as a table in a .mat file.
%   This metadata includes information about the experiment (e.g. time
%   stamps, conditions, camera, temperature), the video (e.g. pixels, frame
%   rate), the behavioral arena (e.g. diameter, mask, orientation), and the
%   lawns (e.g. mask, size, spacing).
%
%   FUNCTIONS:
%
%       GET METADATA
%       - getExperimentInfo: retrieves metadata from a .xlsx file and
%           returns it as info [table]
%       - getVideoInfo: retrieves metadata from a .avi video file and saves
%           it as video [struct]
%
%       GET ARENA PROPERTIES
%       - getAcetateArena: takes an image from a video containing an
%           acetate arena, locates the arena, and calculates information
%           about the video, returned as arena [struct]
%       - getAcetateOrientation: takes an image from a video containing
%           both an acetate arena and a reference dot, locates the
%           reference dot, and calculates information about the orientation
%           of the arena, returned as orientation [struct]
%
%       GET LAWN PROPERTIES
%       - getInvisibleLawns: takes the experimental video and the contrast
%           video (if applicable) and finds/estimates the location of 
%           bacterial lawns in the experiment, returns as lawn [struct]
%       - getManualLawns: updates the lawn locations with any manual
%           corrections that were made in Adobe Photoshop (when needed),
%           returned as lawn [struct]
%
%       CREATE GRIDS
%       - createHexagonalGrid: creates a grid of isometric hexagons and
%           corresponding triangles (6 per hexagon) given information about
%           their size and spacing, returned as hexagon [struct] and
%           triangle [struct]
%       - createDiamondGrid: creates a grid of diamonds and corresponding 
%           triangles (2 per diamond) given information about their size 
%           and spacing, returned as diamond [struct] and triangle [struct]
%       - transformLabeledImage: takes in an image and its properties as a
%           structure and performs given image transformations (i.e.
%           registrations) returned as new images [struct]; performed on the
%           grids (hexagon, diamond, triangle) to register them to arena
%
%   Written 2/8/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEFORAGING, GETEXPERIMENTINFO, GETVIDEOINFO, 
%   GETACETATEARENA, GETACETATEORIENTATION, GETINVISIBLELAWNS,
%   GETMANUALLAWNS, CREATEHEXAGONALGRID, CREATEDIAMONDGRID, 
%   TRANSFORMLABELEDIMAGE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Prepare Workspace and Get Experiments to Analyze

% Add paths
addpath(genpath('Z:\jhaley\foragingPaper\'))

% Set info file (if not specified as an input)
[infoFile,path] = uigetfile({'*.xls*','*.xlsx'});

% Set directory to save all files in
saveDir = [path,'videos',filesep];

% Get experiment name
i = [0,strfind(infoFile,filesep)]; % index where the file name starts
expName = infoFile(i(end)+1:strfind(infoFile,'.')-1);

% Get experiment info from notebook file
info = getExperimentInfo([path,infoFile]);

% Get experiment ids
expToAnalyze = unique([info.expNum]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Get Video, Arena, and Lawn Information

% Query user to overwrite already saved data (if applicable)
overwrite = questdlg('Would you like to overwrite the saved arena and lawn properties?',...
    'Overwrite','Yes','No','No');

arenaDiameter = inputdlg('Diameter of arena?','Arena Diameter',[1 45],{'30'});
arenaDiameter = str2double(arenaDiameter{1});

fprintf('Retrieving Experiment Info: \n')

% v = VideoWriter([path,'lawns'],'MPEG-4');
% open(v);
% figure('Position',[500 10 1000 1000]);
% set(gcf,'Toolbar','none')

for i = 1:length(expToAnalyze)

    tic
    fprintf('\t[Experiment %.2d]: [Plate',expToAnalyze(i))

    % Get plate numbers for this experiment
    plateNums = info.plateNum([info.expNum] == expToAnalyze(i));
    videoNums = info.videoNum([info.expNum] == expToAnalyze(i));
    numPlates = length(plateNums); % # of plates in this experiment
    expDate = datestr(mean(info.timeRecord(info.expNum == expToAnalyze(i)),...
        'omitnan'),'yy-mm-dd');
    
    % Check if arena folder exists to save metadata
    arenaFolder = [saveDir,expDate,'\data\arena\'];
    if ~exist(arenaFolder,'dir')
        mkdir(arenaFolder)
    end
    
    % Check if lawn folder exists to save metadata
    lawnFolder = [saveDir,expDate,'\data\lawn\'];
    if ~exist(lawnFolder,'dir')
        mkdir(lawnFolder)
    end
    
    for j = 1:numPlates
        plateNum = plateNums(j);
        videoNum = videoNums(j);
        ind = info.plateNum == plateNum & info.videoNum == videoNum;
        
        % Build file path to save analyzed arena
        arenaFile = [arenaFolder,info.videoFileName{ind}(1:end-3),'mat'];
        
        % Get video and arena information
        if exist(arenaFile,'file') && strcmp(overwrite,'No')
            warning('off','MATLAB:audiovideo:VideoReader:FileNotFound');
            load(arenaFile,'video','arena','orientation','hexagon','diamond','triangle');
        elseif isempty(info.videoFileName{ind})
            continue
        else
            try
                video = getVideoInfo(info.videoFileName{ind});
            catch ME
                if strcmp(ME.identifier,'MATLAB:audiovideo:VideoReader:UnknownCodec')
                    continue
                end
            end
            arena = getAcetateArena(video.firstFrame,arenaDiameter);
            if contains(path,'Mini')
                orientation.mask = false(size(arena.mask));
                orientation.center = [0,0];
                orientation.distance = NaN;
                orientation.angle = 0;
            else
                orientation = getAcetateOrientation(video.firstFrame,arena);
            end
            [hexagon,~] = createHexagonalGrid(size(arena.mask),...
                arena.scale,arena.diameter,info.lawnDiameter(ind),...
                info.lawnSpacing(ind));
            [diamond,triangle] = createDiamondGrid(size(arena.mask),...
                arena.scale,arena.diameter,info.lawnDiameter(ind),...
                info.lawnSpacing(ind));
            hexagon = transformLabeledImage(hexagon,arena.offset,-orientation.angle,1);
            diamond = transformLabeledImage(diamond,arena.offset,-orientation.angle,1);
            triangle = transformLabeledImage(triangle,arena.offset,-orientation.angle,1);
            save(arenaFile,'video','arena','orientation','hexagon','diamond','triangle');
        end
        
        % Get video object, # of frames, and first frame
        info.numFrames(ind) = video.numFrames;
        info.frameRate(ind) = video.frameRate;
        info.pixels(ind) = {video.pixels};
        info.firstFrame(ind) = {video.firstFrame};
        
        % Get arena and orientation mask, scale, shape, size, and position
        info.scale(ind) = arena.scale;
        info.arenaMask(ind) = {arena.mask};
        info.arenaCenter(ind) = {arena.center};
        info.arenaOffset(ind) = {arena.offset};
        info.arenaAngle(ind) = orientation.angle;
        info.arenaDiameter(ind) = arena.diameter;
        info.arenaShape(ind) = {arena.shape};
        info.arenaCircularity(ind) = arena.circularity;
        info.refMask(ind) = {orientation.mask};
        info.refCenter(ind) = {orientation.center};
        info.refDistance(ind) = orientation.distance;

        % Get grid labels, centers, and area
        info.hexagonGrid(ind) = {hexagon.labeled};
        info.hexagonCenters(ind) = {hexagon.centers};
        info.hexagonArea(ind) = {hexagon.area};
        info.diamondGrid(ind) = {diamond.labeled};
        info.diamondCenters(ind) = {diamond.centers};
        info.diamondArea(ind) = {diamond.area};
        info.triangleGrid(ind) = {triangle.labeled};
        info.triangleCenters(ind) = {triangle.centers};
        info.triangleArea(ind) = {triangle.area};

        % Build file path to save analyzed lawn
        if ~isempty(info.lawnFileName{ind})
            lawnFile = [lawnFolder,info.lawnFileName{ind}(1:end-3),'mat'];
        else
            lawnFile = [lawnFolder,info.videoFileName{ind}(1:end-3),'mat'];
        end
        
        % Get lawn(s) and statistics
        if exist(lawnFile,'file') && strcmp(overwrite,'No')
            load(lawnFile,'lawn');
            % If a manually adjusted file exists, use that
            if exist([lawnFile(1:end-4),'_manual.png'],'file')
                lawn = getManualLawns([lawnFile(1:end-4),'_manual.png'],lawn,arena.scale);
                save(lawnFile,'lawn');
            elseif exist([lawnFile(1:end-4),'_mask.png'],'file')
                lawn = getManualLawns([lawnFile(1:end-4),'_mask.png'],lawn,arena.scale);
                save(lawnFile,'lawn');
            end
        else
            lawn = getInvisibleLawns(info.lawnFileName{ind},...
                video.firstFrame,arena,orientation,...
                info.lawnDiameter(ind),info.lawnSpacing(ind),info.OD600(ind));
            save(lawnFile,'lawn');
        end

        % Add lawn properties to info
        info.lawnCenters(ind) = {lawn.centers};
        info.lawnRadii(ind) = {lawn.radii};
        info.lawnArea(ind) = {lawn.area};
        info.lawnDiameterMean(ind) = {lawn.diameterMean};
        info.lawnDiameterSTD(ind) = {lawn.diameterSTD};
        info.lawnDistributionCenter(ind) = {lawn.distributionCenter};
        info.lawnDistributionSTD(ind) = {lawn.distributionSTD};
        info.lawnSpacingMean(ind) = {lawn.spacingMean};
        info.lawnSpacingSTD(ind) = {lawn.spacingSTD};
        info.lawnCircularity(ind) = {lawn.circularity};
        info.lawnMask(ind) = {lawn.mask};
        info.lawnClosest(ind) = {lawn.closest};
        info.lawnClosestOD600(ind) = {lawn.closestOD600};
        info.lawnMethod(ind) = {lawn.method};
        info.lawnRegistration(ind) = lawn.regFit;

        fprintf(' %.3d',plateNum)
    end
    fprintf('] %.1f s \n',toc)
end
% writeVideo(v,frame);
% close(v)

% Save info file
save([path,expName,'.mat'],'info','-v7.3');
fprintf('Info saved: %s \n',[path,expName,'.mat'])

end