function [info] = getPlateInfo(infoFileName)
% [info] = GETPLATEINFO(infoFileName)
%
%   GETPLATEINFO takes the .xlsx document (infoFileName) containing 
%   metadata information about bacterial plates and returns a table (info).
%   This function is intended for use with GFP bacterial imaging
%   experiments (i.e. plates with no animals).
%
%   INPUTS:
%       - infoFileName [str]: path to a specially formated .xlsx file
%
%   OUTPUTS:
%       - info [table]: a table containing the following variables
%           - plateNum [double]: unique number identifier of each 
%               experimental plate (e.g. 1)
%           - expNum [double]: unique number identifier of each 
%               experimental date; plates on this date share the same
%               growth conditions (e.g. 1)
%           - template [1xN char]: name of the acetate template used (e.g.
%               'rectangle', 'arena', or 'none')
%           - peptone [str]: describes whether the condition plate was made
%               with or without peptone (i.e. 'with' or 'without')
%           - timePoured [datetime]: time that the experimental plates were
%               poured (i.e. DD-Mon-YYYY HH:MM:SS)
%           - timePouredColdRoom [datetime]: time that the experimental 
%               plates were put in the cold room after being poured
%               (i.e. DD-Mon-YYYY HH:MM:SS)
%           - plateAgarDryTimeRoomTemp [duration]: duration of time agar was
%               allowed to dry at room temperature (i.e. HH:MM:SS)
%           - timeSeed [datetime]: time that the experimental plates were
%               seeded (i.e. DD-Mon-YYYY HH:MM:SS)
%           - timeSeedColdRoom [datetime]: time that the experimental plates
%               were put in the cold room after being seeded
%               (i.e. DD-Mon-YYYY HH:MM:SS)
%           - plateBacteriaDryTimeRoomTemp [duration]: duration of time 
%               bacterial lawns were allowed to dry at room temperature 
%               (i.e. HH:MM:SS)
%           - timeRoomTemp [datetime]: time that the experimental plates
%               were removed from the 4C cold room and allowed to warm and
%               grow bacteria (i.e. DD-Mon-YYYY HH:MM:SS)
%           - bacteria [str]: name of bacterial strain (e.g. 'OP50')
%           - OD600Real [double]: measured estimate of the OD600 of the
%               master solution with expected value of 10 (e.g. 10.2875); 
%               used to assess variability of pipetted solutions
%           - CFU [double]: measured value of the CFU count for a dilution
%               series (e.g. 580); used to assess variability of pipetted
%               solutions
%           - OD600 [double]: intended OD600 of solutions pipetted on each 
%               plate (e.g. 1)
%           - lawnVolume [double]: pipetted volume of each lawn in uL 
%               (e.g. 0.5)
%           - lawnDiameter [double]: expected diameter in mm of each lawn 
%               given the volume of liquid pipetted (e.g. 1.5)
%   
%   Written 4/13/2022 by Jess Haley in MATLAB R2021a.
%
%   See also GETEXPERIMENTINFO, ANALYZEGFP, READTABLE, DATETIME, DATESHIFT, 
%   DURATION, DAYS, STR2DOUBLE, TABLE2CELL, CELL2STRUCT.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn off warnings
warning('off','MATLAB:table:RowsAddedExistingVars')

% Define readtable() options and import as a structure
opts = detectImportOptions(infoFileName);
opts.DataRange = 'A1';
opts.VariableTypes = repmat({'char'},1,length(opts.VariableTypes));
experimentInfo = table2cell(...
    readtable(infoFileName,...
    opts,'ReadRowNames',false,'ReadVariableNames',false));
experimentInfo = cell2struct(experimentInfo(:,2:end),...
    experimentInfo(:,1),1);

% Add metadata to a table
info = table();
for i = 1:size(experimentInfo,1)
    % Assign Plate Number
    info.plateNum(i) = str2double(experimentInfo(i).plateNumber);
    
    % Get Experiment Number
    info.expNum(i) = str2double(experimentInfo(i).experimentNumber);

    % Get Template Type
    info.template(i) = {experimentInfo(i).templateType};

    % Get Plate Type
    info.peptone(i) = {experimentInfo(i).platePeptone};

    % Get Date + Time Plate was Poured
    if ~isempty(experimentInfo(i).platePouredDay)
        info.timePoured(i) = dateshift(datetime([experimentInfo(i).platePouredDay,' ',...
            char(duration([24*str2double(experimentInfo(i).platePouredTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    end

    % Get Date + Time Plate was put in Cold Room after Drying
    if ~isempty(experimentInfo(i).platePouredColdRoomDay)
        info.timePouredColdRoom(i) = dateshift(datetime([experimentInfo(i).platePouredColdRoomDay,' ',...
            char(duration([24*str2double(experimentInfo(i).platePouredColdRoomTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
        info.plateAgarDryTimeRoomTemp(i) = info.timePouredColdRoom(i)-info.timePoured(i);
    end
    
    % Get Date + Time Plate was Seeded; Get Duration of Lawn Growth
    if ~isempty(experimentInfo(i).plateSeededDay)
        info.timeSeed(i) = dateshift(datetime([experimentInfo(i).plateSeededDay,' ',...
            char(duration([24*str2double(experimentInfo(i).plateSeededTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    end

    % Get Date + Time Plate was put in Cold Room after Seeding
    if ~isempty(experimentInfo(i).plateSeededColdRoomDay)
        info.timeSeedColdRoom(i) = dateshift(datetime([experimentInfo(i).plateSeededColdRoomDay,' ',...
            char(duration([24*str2double(experimentInfo(i).plateSeededColdRoomTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
        info.plateBacteriaDryTimeRoomTemp(i) = info.timeSeedColdRoom(i)-info.timeSeed(i);
    end
    
    % Get Date + Time Plate was Removed from Cold Room for Imaging
    if ~isempty(experimentInfo(i).plateRoomTempTime)
        info.timeRoomTemp(i) = dateshift(datetime([experimentInfo(i).plateRoomTempDay,' ',...
            char(duration([24*str2double(experimentInfo(i).plateRoomTempTime),0,0]))],...
            'InputFormat','dd-MMM-yyyy HH:mm:ss'),'start','minute','nearest');
    end
    
    % Get Bacterial Identifier & Estimated CFUs (in 0.5 uL of OD600 = 10)
    info.bacteria(i) = {experimentInfo(i).bacteriaSpecies};
    info.OD600Real(i) = str2double(experimentInfo(i).bacteriaOD600);
    info.CFU(i) = str2double(experimentInfo(i).bacteriaCFU);
    
    % Get Bacterial Concentration Seeded (OD600)
    info.OD600(i) = str2double(experimentInfo(i).plateSeededOD600);

    % Get Volume (uL) and Intendended Diameter (mm) per Lawn
    info.lawnVolume(i) = str2double(experimentInfo(i).plateSeededVolume);
    info.lawnDiameter(i) = str2double(experimentInfo(i).plateSeededLawnDiameter);
    
end

end