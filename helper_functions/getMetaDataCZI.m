function [metaData,images] = getMetaDataCZI(folder)
% [metaData] = GETMETADATACZI(folder)
%
%   GETMETADATACZI takes a directory (folder) containing Zeiss image files
%   (.czi) and saves the metadata as a table (metadata.csv) and the images
%   as a cell array (data.mat) in the directory.
%
%   INPUTS:
%       - folder [str]: path to a directory containing .czi files (or a
%           directory containing subfolders containing .czi files)
%               
%   OUTPUTS:
%       - metaData [table]: a table containing the following variables
%           - fileName [str]: file name of the Zeiss image file (*.czi)
%           - pathName [str]: file path where the image file is saved
%           - savedTime [datetime]: time that the image was saved 
%               (i.e. DD-Mon-YYYY HH:MM:SS); time zone is 'local'
%           - acquisitionTime [datetime]: time that the image was captured 
%               (i.e. DD-Mon-YYYY HH:MM:SS); time zone is 'local'
%           - exposureTime [double]: image exposure time in ms (e.g. 300)
%           - bitDepth [str]: image data type (e.g. uint16)
%           - xPixels [double]: image width in pixels (e.g. 1200)
%           - yPixels [double]: image height in pixels (e.g. 1200)
%           - minValue [num]: minimum pixel intensity of image
%           - maxValue [num]: maximum pixel intensity of image
%           - meanValue [num]: mean pixel intensity of image
%           - fluorescence [logical]: estimate of if the image was
%               fluorescence (1) or brightfield (0) microscopy
%       - images [cell]: a cell array containing each image from the input
%           directory
%   
%   Written 1/24/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEGFP, BFOPEN, DIR, DATETIME, CELLFUN, CLASS.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn off warnings
warning('off','MATLAB:table:RowsAddedExistingVars')

% Get files
files = dir([folder,filesep,'*.czi*']);
if isempty(files)
    files = dir([folder,filesep,'**',filesep,'*.czi*']);
end
numFiles = size(files,1);

% Get file names and paths
fileName = {files(:).name}';
pathName = {files(:).folder}';

% Get the date and time that the image was saved
savedTime = cellfun(@(d) datetime(d,'TimeZone','local'),{files(:).date}');

% Initialize variables
images = cell(numFiles,1);
acquisitionTime = NaT(numFiles,1,'TimeZone','UTC');
exposureTime = nan(numFiles,1);

% Start parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

parfor i = 1:numFiles

    % Use bioformats package to open the file
    bfData = bfopen([files(i).folder,filesep,files(i).name]);

    % Get the image and metadata info
    images{i} = bfData{1,1}{1};
    metadata = bfData{1,2};
    omeMeta = bfData{1,4};

    % Get the date and time that the image was acquired (in local time)
    acquisitionTime(i) = datetime(char(omeMeta.getImageAcquisitionDate(0)),...
        'InputFormat','uuuu-MM-dd''T''HH:mm:ss','TimeZone','UTC');

    % Get the exposure time (ms)
    exposureTime(i) = str2double(metadata.get('Global HardwareSetting|ParameterCollection|ExposureTime'));
end

% Change acquisition time zone
acquisitionTime.TimeZone = 'local';

% Get image info (i.e. bitDepth, height, width)
bitDepth = cellfun(@(i) class(i),images,'UniformOutput',false);
xPixels = cellfun(@(i) size(i,1), images);
yPixels = cellfun(@(i) size(i,2), images);

% Get mean, max, and min pixel intensity
minValue = cellfun(@(i) min(i,[],'all'),images);
maxValue = cellfun(@(i) max(i,[],'all'),images);
meanValue = cellfun(@(i) mean(i,'all'),images);

% Determine if images are brightfield (0) or fluorescence (1)
fracMax = meanValue./double((cellfun(@(i) intmax(i),bitDepth))); % fraction of maximum pixel value
fluorescence = fracMax./exposureTime < 1e-3;

% Add info to the table
metaData = table(fileName,pathName,savedTime,acquisitionTime,...
    exposureTime,bitDepth,xPixels,yPixels,minValue,maxValue,...
    meanValue,fluorescence);

% Write the table to file
writetable(metaData,[folder,filesep,'metadata.csv'])
save([folder,filesep,'data.mat'],'images','-v7.3')

end