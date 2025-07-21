function [info] = getExperimentInfo(infoFileName)
% This function takes a .xlsx document containing experimental metadata
% information and creates an info table.
%   INPUTS:
%       - infoFileName {str}: path to a specially formated .xlsx file
%   OUTPUTS:
%       - info [table]: a table containing the following variables
%           - plateNum [uint16]: unique number identifier of each 
%               experimental plate; worms on these plates share the same
%               growth and experimental conditions (e.g. 1)
%           - videoNum [uint16]: unique number identifier for each video
%               recorded for each plate (e.g. 1); only applicable for
%               experiments where multiple videos were recorded
%           - expNum [uint16]: unique number identifier of each 
%               experimental date; plates on this date share the same
%               growth conditions (e.g. 1)
%           - wormNum {[1xN uint16]}: unique number identifier(s) for all
%               worm(s) on each experimental plate (e.g. [1,2,3,4])
%           - condition {1xN char}: name of the condition tested (e.g.
%               'grid' or 'single')
%           - strainName [str]: unique identifier for the name of the
%               strain (e.g. well-fed, food-deprived, or osm-6)
%           - strainID [str]: unique identifier for the CGC strain id (e.g.
%               N2, or PR811)
%           - exclude [logical]: true if experiments are to be excluded
%               from future analyses (e.g. experiments lacking a contrast
%               video where lawns cannot be reliably detected)
%           - timeRecord [datetime]: time that the experiemntal recording 
%               started (i.e. DD-Mon-YYYY HH:MM:SS)
%           - growthTimePicked [datetime]: time that L4 animals were picked
%               onto a "growth" plate (where they develop from L4 -> 1DOA)
%               on the day(s) preceeding the experiment (i.e. 
%               DD-Mon-YYYY HH:MM:SS)
%           - growthTimeSeed [datetime]: time that the growth plates were
%               seeded (i.e. DD-Mon-YYYY HH:MM:SS)
%           - growthTimeColdRoom [datetime]: time that the growth plates
%               were placed in the 4C cold room after seeding
%               (i.e. DD-Mon-YYYY HH:MM:SS)
%           - growthTimeRoomTemp [datetime]: time that the growth plates
%               were removed from the 4C cold room and allowed to warm and
%               grow bacteria (i.e. DD-Mon-YYYY HH:MM:SS)
%           - growthlawnGrowth [duration]: duration of time growth plate
%               lawns were allowed to grow for at room temperature after
%               drying (i.e. HH:MM:SS)
%           - growthOD600 [double]: OD600 of dilluted solutions pipetted
%               onto each growth plate (e.g. 1)
%           - peptone [str]: describes whether the condition plate was made
%               with or without peptone (i.e. 'with' or 'without')
%           - timeSeed [datetime]: time that the experimental plates were
%               seeded (i.e. DD-Mon-YYYY HH:MM:SS)
%           - timeColdRoom [datetime]: time that the experimental plates
%               were placed in the 4C cold room after seeding
%               (i.e. DD-Mon-YYYY HH:MM:SS)
%           - timeRoomTemp [datetime]: time that the experimental plates
%               were removed from the 4C cold room and allowed to warm and
%               grow bacteria (i.e. DD-Mon-YYYY HH:MM:SS)
%           - lawnGrowth [duration]: duration of time lawns were allowed to
%               grow for at room temperature after drying (i.e. HH:MM:SS)
%           - wormGrowth [duration]: duration of time between L4s being
%               picked and the experiment start time (i.e. HH:MM:SS)
%           - age [str]: approximate life stage of the worm given
%               wormGrowth (i.e. 'L4', '1DOA', or '2DOA')
%           - bacteria [str]: name of bacterial strain (e.g. 'OP50')
%           - OD600 [double or 1xnLawn array]: intended OD600 of solutions
%               pipetted on each plate (e.g. 1); for experiments with
%               multiple patch concentrations, the OD600 is given as a cell
%               containing a double array with OD600 of each patch (e.g.
%               {[1,5,1,5,1,5,1,5]})
%           - OD600Label [str]: OD600 solution pipetted on each plate (e.g.
%               '1.00' or '0.50'); for experiments with multiple patch 
%               concentrations, solutions are listed (e.g. '1.00 5.00
%               10.00')
%           - OD600Real [double]: measured estimate of the OD600 of the
%               master solution with expected value of 10 (e.g. 10.2875); 
%               used to assess variability of pipetted solutions
%           - CFU [uint16]: measured value of the CFU count for a dilution
%               series (e.g. 580); used to assess variability of pipetted
%               solutions
%           - configuration [str]: for experiments with multiple OD600
%               solutions per plate, the pattern of seeding is designated
%               by a configuration code (i.e. 'A', 'B', or 'C')
%           - growthCondition [double]: hours of growth (i.e. 0, 12, or 48)
%           - lawnVolume [double]: pipetted volume of each lawn in uL 
%               (e.g. 0.5)
%           - lawnDiameter [double]: expected diameter in mm of each lawn 
%               given the volume of liquid pipetted (e.g. 1.5)
%           - lawnSpacing [double]: intended center-to-center spacing in mm
%               of each isometric lawn (e.g. 6)
%           - videoFileName {str}: file name of the .avi experimental 
%               recording
%           - lawnFileName {str}: file name of the .avi contrast video 
%               taken to estimate locations of dilute lawns (if available)
%           - wormLabFileName {str}: file name of the .csv exported with
%               the mid-point locations of worms in each plate after
%               tracking in WormLab
%           - temp [double]: average temperature of the room during
%               recording of the experiment in C (e.g. 21.68)
%           - humidity [double]: average percent humidity of the room 
%               during recording of the experiment (e.g. 47.71)
%           - camera [uint16]: camera that each recording was taken on
%               (e.g. 1 or 2); each camera can be assumed to have the same
%               field-of-view and focus on a given experimental day
%   
%   Written 2/8/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEFORAGING, READMATRIX, STR2DOUBLE, DATETIME, DATESHIFT,
%   DURATION, DAYS, MINUTE, CHAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn off warnings
warning('off','MATLAB:table:RowsAddedExistingVars')

% Define readtable() options and import as a structure
opts = detectImportOptions(infoFileName);
opts.DataRange = 'A1';
opts.VariableTypes = repmat({'char'},1,length(opts.VariableTypes));
experimentInfo = table2cell(readtable(infoFileName,...
    opts,'ReadRowNames',false,'ReadVariableNames',false));
experimentInfo = cell2struct(experimentInfo(:,2:end),...
    experimentInfo(:,1),1);

% Add plate numbers to structure
plateNum = num2cell(1:size(experimentInfo,1));
[experimentInfo.plateNum] = plateNum{:};
[experimentInfo.videoNum] = deal(1);

% Check for 2nd video file and add rows if necessary
if isfield(experimentInfo,'conditionVideo2')
    experimentInfo2 = experimentInfo;
    [experimentInfo2.conditionVideo] = experimentInfo2(:).conditionVideo2;
    include = false(size(experimentInfo2,1),1);
    for i = 1:size(experimentInfo2,1)
        include(i) = ~isempty(experimentInfo2(i).conditionVideo);
    end
    [experimentInfo2.videoNum] = deal(2);
    experimentInfo = struct2table([experimentInfo;experimentInfo2(include,:)]);
    experimentInfo = table2struct(sortrows(experimentInfo,'plateNum'));
end

% Add metadata to a table
info = table();
for i = 1:size(experimentInfo,1)
    % Assign Plate Number
    info.plateNum(i) = uint16(experimentInfo(i).plateNum);

    % Assign Plate Number
    info.videoNum(i) = uint16(experimentInfo(i).videoNum);
    
    % Get Experiment Number
    info.expNum(i) = uint16(str2double(experimentInfo(i).experimentNumber));
    
    % Get Worm Identification Number(s)
    if length(experimentInfo(i).wormNumber) == 4
        info.wormNum(i) = {str2double(experimentInfo(i).wormNumber)};
    else
        info.wormNum(i) = {str2double(experimentInfo(i).wormNumber(1:4)):...
            str2double(experimentInfo(i).wormNumber(6:end))};
    end
    
    % Get Condition Name
    info.condition(i) = {experimentInfo(i).conditionName};

    % Get Strain Name + ID
    if isfield(experimentInfo,'strainName')
        info.strainName(i) = {experimentInfo(i).strainName};
        info.strainID(i) = {experimentInfo(i).strainID};
    else
        info.strainName(i) = {'well-fed'};
        info.strainID(i) = {'N2'};
    end

    % Get Exclude
    if isfield(experimentInfo,'exclude')
        info.exclude(i) = strcmp(experimentInfo(i).exclude,'exclude');
    else
        info.exclude(i) = false;
    end
    
    % Get Date + Time of Video
    info.timeRecord(i) = datetime(experimentInfo(i).conditionVideo(1:end-6),...
        'InputFormat','yyyy-MM-dd_HH-mm-ss');

    % Get Time when L4 was Picked
    if ~isempty(experimentInfo(i).pickedL4Day)
        info.growthTimePicked(i) = dateshift(datetime([experimentInfo(i).pickedL4Day,' ',...
            char(duration([24*str2double(experimentInfo(i).pickedL4Time),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    else
        info.growthTimePicked(i) = NaT;
    end

    % Get Date + Time L4 Plate was Seeded; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).pickedL4PlateSeededTime)
        info.growthTimeSeed(i) = dateshift(datetime([experimentInfo(i).pickedL4PlateSeededDay,' ',...
            char(duration([24*str2double(experimentInfo(i).pickedL4PlateSeededTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    else
        info.growthTimeSeed(i) = NaT;
    end

    % Get Date + Time L4 Plate was Placed at 4C; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).pickedL4PlateColdRoomTime)
        info.growthTimeColdRoom(i) = dateshift(datetime([experimentInfo(i).pickedL4PlateColdRoomDay,' ',...
            char(duration([24*str2double(experimentInfo(i).pickedL4PlateColdRoomTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    else
        info.growthTimeColdRoom(i) = info.growthTimeSeed(i) + minutes(10);
    end

    % Get Date + Time Growth Plate was Placed at Room Temp; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).pickedL4PlateRoomTempTime)
        info.growthTimeRoomTemp(i) = dateshift(datetime([experimentInfo(i).pickedL4PlateRoomTempDay,' ',...
            char(duration([24*str2double(experimentInfo(i).pickedL4PlateRoomTempTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    else
        info.growthTimeRoomTemp(i) = info.growthTimePicked(i) - minutes(60);
    end
    info.growthLawnGrowth(i) =  info.growthTimePicked(i) - info.growthTimeRoomTemp(i) + ...
        info.growthTimeColdRoom(i) - info.growthTimeSeed(i);

    % Get Growth OD600
    info.growthOD600(i) = str2double(experimentInfo(i).pickedL4PlateOD600);

    % Get Starvation Info (if applicable)
    if isfield(experimentInfo,'starvationDay')
        if ~isempty(experimentInfo(i).starvationDay)
            % Get Time when Worm was Picked to Starvation Plate
            info.starvedTime(i) = dateshift(datetime([experimentInfo(i).starvationDay,' ',...
                char(duration([24*str2double(experimentInfo(i).starvationTime),0,0]))],...
                'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
            
            % Get Duration of Food Deprivation
            info.starvedDuration(i) = info.timeRecord(i) - info.starvedTime(i);
        else
            info.starvedTime(i) = NaT;
            info.starvedDuration(i) = duration(0,0,0);
        end
    end

    % Get Peptone Info
    if isfield(experimentInfo,'conditionPlatePeptone')
        info.peptone(i) = {experimentInfo(i).conditionPlatePeptone};
    else
        info.peptone(i) = {'with'};
    end
    
    % Get Date + Time Plate was Seeded; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).conditionPlateSeededTime)
        info.timeSeed(i) = dateshift(datetime([experimentInfo(i).conditionPlateSeededDay,' ',...
            char(duration([24*str2double(experimentInfo(i).conditionPlateSeededTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    end

    % Get Date + Time Plate was Placed at 4C; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).conditionPlateColdRoomTime)
        info.timeColdRoom(i) = dateshift(datetime([experimentInfo(i).conditionPlateColdRoomDay,' ',...
            char(duration([24*str2double(experimentInfo(i).conditionPlateColdRoomTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    else
        info.timeColdRoom(i) = info.timeSeed(i) + minutes(10);
    end
    
    % Get Date + Time Plate was Placed at Room Temp; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).conditionPlateRoomTempTime)
        info.timeRoomTemp(i) = dateshift(datetime([experimentInfo(i).conditionPlateRoomTempDay,' ',...
            char(duration([24*str2double(experimentInfo(i).conditionPlateRoomTempTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    else
        info.timeRoomTemp(i) = info.timeRecord(i) - minutes(60);
    end
    info.lawnGrowth(i) = info.timeRecord(i)-info.timeRoomTemp(i) + ...
        info.timeColdRoom(i) - info.timeSeed(i);
    
    % Get Worm Growth
    if ~isempty(info.growthTimePicked(i))
        info.wormGrowth(i) = info.timeRecord(i)-info.growthTimePicked(i);
    end

    % Get Worm Age
    ageInDays = days(datetime(experimentInfo(i).experimentDay,'InputFormat','dd-MMM-yyyy') - ...
        datetime(experimentInfo(i).pickedL4Day,'InputFormat','dd-MMM-yyyy'));
    ageStrings = {'L4','1DOA','2DOA'};
    info.age(i) = ageStrings(ageInDays + 1);
    
    % Get Bacterial Identifier & Estimated CFUs (in 0.5 uL of OD600 = 10)
    info.bacteria(i) = {experimentInfo(i).conditionPlateBacteria};
    info.OD600Real(i) = str2double(experimentInfo(i).conditionPlateOD600Real);
    info.CFU(i) = str2double(experimentInfo(i).conditionPlateCFU);
    
    % Get Bacterial Concentration (OD600)
    if isfield(experimentInfo,'conditionPlateOD600')
        slash = strfind(experimentInfo(i).conditionPlateOD600,'/');
        if isempty(experimentInfo(i).conditionPlateOD600)
            info.OD600(i) = 0;
        elseif isempty(slash)
            info.OD600(i) = str2double(experimentInfo(i).conditionPlateOD600);
        else
            info.OD600(i) = [str2double(experimentInfo(i).conditionPlateOD600(1:slash-1)),...
                str2double(experimentInfo(i).conditionPlateOD600(slash+1:end))];
        end
        info.OD600Label(i) = {num2str(info.OD600(i),'% 0.2f')};
    else
        info.configuration(i) = {experimentInfo(i).conditionConfiguration};
        OD = [str2double(experimentInfo(i).conditionPlateOD600_1),...
            str2double(experimentInfo(i).conditionPlateOD600_2),...
            str2double(experimentInfo(i).conditionPlateOD600_3)];
        switch info.configuration{i}
            case 'A'
                info.OD600(i) = {[OD([2,1,1,2,2,2,1,1,1]),0,OD([1,2,2,2,1,2,1,1,2])]};
            case 'B'
                info.OD600(i) = {[OD([1,2,2,1,2,1,1,1,2]),0,OD([2,2,2,1,1,1,2,2,1])]};
            case 'C'
                info.OD600(i) = {[OD([2,3,1,1,3,3,2,1,2]),0,OD([2,1,2,3,3,1,1,3,2])]};
        end
        info.OD600Label(i) = {num2str(unique(OD(~isnan(OD))),'% 0.2f')};
    end

    % Get Growth Condition
    info.growthCondition(i) = 12*round(hours(info.lawnGrowth(i))/12);
    
    % Get Volume (uL) per Lawn
    if isempty(experimentInfo(i).conditionPlateVolume)
        info.lawnVolume(i) = 0;
    else
        info.lawnVolume(i) = str2double(experimentInfo(i).conditionPlateVolume);
    end
    
    % Get Intended Lawn Diameter (mm)
    if isempty(experimentInfo(i).conditionPlateLawnDiameter)
        info.lawnDiameter(i) = 0;
    else
        info.lawnDiameter(i) = str2double(experimentInfo(i).conditionPlateLawnDiameter);
    end
    
    % Get Lawn Center-to-Center Spacing (mm)
    if isempty(experimentInfo(i).conditionPlateLawnSpacingMean)
        info.lawnSpacing(i) = 0;
    else
        info.lawnSpacing(i) = str2double(experimentInfo(i).conditionPlateLawnSpacingMean);
    end
    
    % Get File Names
    info.videoFileName(i) = {experimentInfo(i).conditionVideo};
    info.lawnFileName(i) = {experimentInfo(i).conditionContrastVideo};
    if ~strcmp(info.videoFileName(i),'')
        info.wormLabFileName(i) = {[experimentInfo(i).conditionVideo(1:end-4),'.csv']};
    else
        info.wormLabFileName(i) = {''};
    end
    
    % Get Temperature (C) + Humidity (%)
    if ~isempty(experimentInfo(i).conditionTemp)
        info.temp(i) = str2double(experimentInfo(i).conditionTemp);
        info.humidity(i) = str2double(experimentInfo(i).conditionHumidity);
    else
        info.temp(i) = NaN;
        info.humidity(i) = NaN;
    end

    % Get Experiment Camera Number (1 or 2)
    if ~isempty(experimentInfo(i).conditionVideo)
        info.camera(i) = str2double(experimentInfo(i).conditionVideo(end-4));
    end
end

end