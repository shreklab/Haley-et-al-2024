function [] = createVideo(data,info,saveDir)
%[] = CREATEVIDEO(data,info)
%
%   CREATEVIDEO writes a video to file showing either the tracks of an
%   animal or the original video downsampled with scale bar and time stamp
%   showing.
%
%   INPUTS:
%       - data [struct]: structure containing the location of worm(s) at 
%           every timepoint; created with ANALYZEFORAGING
%       - info [struct]: structure containing the metadata for all worm(s);
%           created with GETFORAGINGINFO
%       - saveDir [str]: path name to save the video file
%
%   Written 5/8/2023 by Jess Haley in MATLAB R2023a.
%
%   See also ANALYZEFORAGING, GETFORAGINGINFO, PLOTTRACKS.

% Set properties
prompts = {'wormNum:','length of video (seconds):',...
    'video content ("track" or "original"):','arena mask ("yes" or "no"):',...
    'scale bar size (mm):','simulate eggs ("yes" or "no"):'};
defaults = {'1','60','track','yes','5','no'};
answers = inputdlg(prompts,'Video Properties',[1 50],defaults);
wormNum = str2double(answers{1});
videoLength = str2double(answers{2}); % seconds
content = answers{3}; % 'track' or 'original'
mask = answers{4}; % 'yes' or 'no'
scaleBarSize = str2double(answers{5}); % mm
addEggs = answers{6}; % simulation to add fake eggs

% Get plate and video number
plateNum = data.plateNum(find(data.wormNum == wormNum,1));
videoNum = data.videoNum(find(data.wormNum == wormNum,1));
expNum = data.expNum(find(data.wormNum == wormNum,1));
indInfo = info.plateNum == plateNum & info.videoNum == videoNum;
indData = data.wormNum == wormNum;

% Get border amplitudes to use for patch color
borderAmplitude = readtable('Z:\jhaley\foragingPaper\foragingGFP\borderAmplitude.csv');
uniqueOD600 = unique(info.lawnClosestOD600{indInfo});
if uniqueOD600 > 1e-2
    indBorder = find(strcmp(borderAmplitude.peptone,info.peptone(indInfo)) & ...
        borderAmplitude.growthTimeCondition >= info.growthCondition(indInfo) & ...
        ismember(borderAmplitude.OD600,uniqueOD600) & borderAmplitude.lawnVolume == info.lawnVolume(indInfo));
    ind10 = find(strcmp(borderAmplitude.peptone,'with') & ...
        borderAmplitude.growthTimeCondition == 0 & ...
        borderAmplitude.OD600 == 10 & borderAmplitude.lawnVolume == 0.5);
    borderValue = borderAmplitude.slope(indBorder).*60.*(info.growthCondition(indInfo)+1) + ...
        borderAmplitude.intercept(indBorder);
    borderValue10 = borderAmplitude.slope(ind10).*92 + borderAmplitude.intercept(ind10);
    borderValue = 10*borderValue./borderValue10;
else
    borderValue = 1e-2;
end
colorValue = min(1 - (log(borderValue)*0.09 + 0.4),0.8);

% Set video file properties
if nargin >= 3
    videoFolder = saveDir;
else
    videoFolder = [fileparts(which(info.videoFileName{indInfo})),...
        filesep,'plot',filesep,content,'Video',filesep];
    if ~exist(videoFolder,'dir')
        mkdir(videoFolder)
    end
end
videoName = [videoFolder,...
    'Exp',num2str(expNum,'%.4d'),'_Plate',num2str(plateNum,'%.4d'),...
     '_Worm',num2str(wormNum,'%.4d')];
v = VideoWriter(videoName,'MPEG-4');
open(v);

% Get data
wormData = data(indData,:);
wormData.xPosition = round(wormData.xPosition);
wormData.yPosition = round(wormData.yPosition);
frameRate = round(info.frameRate(indInfo));
pixels = info.pixels{indInfo};
scale = info.scale(indInfo);
video = getVideoInfo(info.videoFileName{indInfo});

% Set track properties
trailFrames = 30*frameRate; % # frames that have a darker colored trail
colors = [linspace(230,128,trailFrames-10),zeros(1,10)];
numSteps = 30*videoLength;
edges = round(linspace(1,max(wormData.frameNum),numSteps));

% Set image properties
lawnOD600 = info.lawnClosestOD600{indInfo}.*info.lawnMask{indInfo};
lawnBorder = ones(size(lawnOD600));
for i = 1:length(uniqueOD600)
    lawnBorder(lawnOD600 == uniqueOD600(i)) = colorValue(i);
end
lawnImage = uint8(255.*lawnBorder);
arenaMask = uint8(info.arenaMask{indInfo});

% Set plot properties
figure('Position',[0 0 560 560])
displayImage = flip(lawnImage);
if strcmp(mask,'yes')
    displayImage = flip(arenaMask).*displayImage;
    displayImage(displayImage == 0) = 32; % Keynote gray (actually 34; but MPEG-4 makes that 38)
end
eggImage = uint8(ones(size(displayImage)));
h = imshow(displayImage);
hold on
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')

% Plot scale bar
scaleBar = plot(24.25*pixels/25 - [0 scaleBarSize*scale],0.75*pixels/25 + [0 0],...
    'Color','k','LineWidth',3);
scaleText = text(24.25*pixels(1)/25 - (scaleBarSize/2)*scale,pixels(2)/25,[num2str(scaleBarSize),' mm'],...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Color','k','FontName','Arial','FontSize',24);

% Plot time (s)
timeText = text(24.25*pixels(1)/25,24.25*pixels(2)/25,datestr(seconds(0),'HH:MM:SS'),...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'Color','k','FontName','Arial','FontSize',24);

if strcmp(mask,'yes')
    % Change text color to white
    set(timeText,'Color',[1 1 1])
    set(scaleBar,'Color',[1 1 1])
    set(scaleText,'Color',[1 1 1])
end

if strcmp(content,'track')

    % Plot lawn diameter (if applicable)
    % lawnText = text(0.75*pixels(1)/25,24.25*pixels(2)/25,...
    %     ['Ø ',num2str(info.lawnDiameterMean{indInfo},'%.1f'),' mm'],...
    %     'HorizontalAlignment','left','VerticalAlignment','top',...
    %     'FontName','Arial','FontSize',24,'Color',scaleText.Color);
    % 
    % % Plot lawn OD600 (if applicable)
    % try
    % OP50Text = text(24.25*pixels(1)/25,24.25*pixels(2)/25,...
    %     ['OD_6_0_0 = ',num2str(info.OD600(indInfo),'%.0f')],...
    %     'HorizontalAlignment','right','VerticalAlignment','top',...
    %     'FontName','Arial','FontSize',24,'Color',scaleText.Color);
    % end

    % Set track filter
    filterRadius = 2;
    structure = strel('disk',filterRadius); structure = structure.Neighborhood;
    [filterX,filterY] = ind2sub(size(structure),find(structure));
    filterX = filterX - filterRadius; filterY = filterY - filterRadius;
end

for j = 1:numSteps
    currentFrame = edges(j);

    switch content
        case 'track'
            currentTrail = min(currentFrame,trailFrames);
            for f = currentFrame-currentTrail+1:currentFrame
                ind = find(round(wormData.timeOffset*...
                    info.frameRate(indInfo)) == f);
                if ~isempty(ind)
                    if ~isnan(wormData.xPosition(ind))
                        xPosition = wormData.yPosition(ind)+filterY;
                        yPosition = wormData.xPosition(ind)+filterX;
                        position = sub2ind(pixels,xPosition,yPosition);
                        displayImage(position) = ...
                            colors(trailFrames - currentFrame + f);
                    end
                end
            end

        case 'original'
            displayImage = flip(read(video.vidObject,currentFrame));
            if strcmp(mask,'yes')
                displayImage = flip(arenaMask).*displayImage;
                displayImage(displayImage == 0) = 32; % Keynote gray (actually 34; but MPEG-4 makes that 38)  
            end
    end
    
    % Update image and time display
    set(h,'CData',displayImage);
    set(timeText,'String',...
        datestr(seconds(currentFrame/frameRate),'HH:MM:SS'));

    % Add fake eggs
    if strcmp(addEggs,'yes')
        ind = find(wormData.frameNum == currentFrame);
        try
            eggImage(wormData.yPosition(ind),wormData.xPosition(ind)) = rand(1) >= 0.1;
        end
        set(h,'CData',displayImage.*imerode(eggImage,strel('disk',3)));
    end
    
    % Save frame
    frame = getframe(gca);
    writeVideo(v,frame);
end

% Close video
close(v)

end