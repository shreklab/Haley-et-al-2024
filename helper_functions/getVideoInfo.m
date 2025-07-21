function [video] = getVideoInfo(fileName)
% [video] = GETVIDEOINFO(fileName)
%
%   GETVIDEOINFO takes in a video file (fileName) and uses VideoReader to
%   retrieve metadata returned as a structure (video)
%
%   INPUTS:
%       - fileName [str]: file name of a video file (e.g. 
%           YYYY-MM-DD_HH-MM-SS_C.avi); C here is used to denote the camera
%           number
%
%   OUTPUTS:
%       - video [struct]: a structure containing the following fields
%           - numFrames [double]: number of frames in the video (e.g.
%               10797)
%           - frameRate [double]: average frames acquired per second (e.g. 
%               2.9991)
%           - duration [double]: number of seconds in the video (e.g. 
%               3599.8) 
%           - pixels [1x2 double]: video width, height in pixels (e.g.
%               [1024,1024])
%           - format [str]: pixel format (e.g. 'RGB24', 'Mono8',
%               'Grayscale')
%           - bitDepth [double]: bits per pixel (e.g. 8, 16)
%           - megaBytes [double]: video size in MB
%           - camera [double]: camera identifier for which the recording
%               was taken (e.g. 1 or 2)
%           - vidObject [1x1 VideoReader]: video object used with read() to
%               get frames
%           - firstFrame [MxN num]: first frame of the video (e.g.
%               [1024x1024 uint8])
%
%   Written 4/13/2022 by Jess Haley in MATLAB R2022a.
%
%   See also ANALYZEFORAGING, VIDEOREADER, READ, STR2DOUBLE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get video object
vidObject = VideoReader(fileName);

% Get video metadata
video.numFrames = vidObject.NumFrames;
video.frameRate = vidObject.FrameRate; % e.g. 30
video.duration = vidObject.Duration;
video.pixels = [vidObject.Width vidObject.Height]; % e.g. [512,512]
video.format = vidObject.VideoFormat; % e.g. 'RGB24', 'Mono8', 'Grayscale'
video.bitDepth = vidObject.BitsPerPixel; % e.g. 8, 16
video.megaBytes = video.bitDepth * video.pixels(1) * video.pixels(2) * ...
    video.numFrames / 8 / 2^20;
video.camera = str2double(fileName(strfind(fileName,'.')-1));
video.vidObject = vidObject;

% Read and store first frame
video.firstFrame = read(vidObject,1);

end