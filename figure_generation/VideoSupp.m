%% Load Data

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';

load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');

saveDir = [path,'figures\VideoSupp\'];

%% Video S1 - Raw Behavior Video

% wormNum = 360
createVideo(data,info,saveDir)

%% Video S2 - Raw Contrast Video

load([path,'foragingConcentration\videos\22-03-18\data\lawn\2022-03-18_12-55-44_2.mat']);
relevantFrames = 79:138;
numFrames = length(relevantFrames);

v = VideoWriter([saveDir,'VideoS2'],'MPEG-4');
v.FrameRate = 15;
open(v);
figure('Position',[0 0 560 560])
displayImage = flip(read(lawn.video.vidObject,relevantFrames(1)));
h = imshow(displayImage);
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')

for f = 1:numFrames
    displayImage = flip(read(lawn.video.vidObject,relevantFrames(f)));
    set(h,'CData',displayImage);

    % Save frame
    frame = getframe(gca);
    writeVideo(v,frame);
    writeVideo(v,frame);
end

% Close video
close(v)

%% Video S3 - Tracking Video

% wormNum = 360
createVideo(data,info,saveDir)

%% Video S4 - Encounter Detection

wormNum = 214;
frameNum = 4366;
load(fullfile(path,'foragingMini','experimentInfo.mat'),'info');
load(fullfile(path,'foragingMini','encounter.mat'),'encounter');
indInfo = find(cellfun(@(w) any(w == wormNum),info.wormNum));
indEncounter = find(strcmp(encounter.expName,'foragingMini') & encounter.wormNum == wormNum);
video = getVideoInfo(info.videoFileName{indInfo});
centerpointsFile = [path,'foragingMini\videos\23-08-25\wormlab\centerpoints\2023-08-25_11-15-33_2.csv'];
wormTrack = readmatrix(centerpointsFile,'Range',[6,1]);
wormTrack = wormTrack(1:4:end,:);
wormX = wormTrack(:,3:2:end);
wormY = wormTrack(:,4:2:end);

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
lawnOD600 = info.lawnClosestOD600{indInfo}.*info.lawnMask{indInfo};
lawnBorder = ones(size(lawnOD600));
for i = 1:length(uniqueOD600)
    lawnBorder(lawnOD600 == uniqueOD600(i)) = colorValue(i);
end
lawnImage = uint8(255.*lawnBorder);
arenaMask = uint8(info.arenaMask{indInfo});

v = VideoWriter([saveDir,'VideoS4'],'MPEG-4');
open(v);

figure('Position',[0 0 560 560])
displayImage = flip(arenaMask.*lawnImage);
displayImage(displayImage == 0) = 34; % Keynote gray
h = imshow(displayImage);
hold on
scaleBarSize = 2; % mm
scale = info.scale(indInfo);
pixels = info.pixels{indInfo};

% Plot scale bar
scaleBar = plot(24.25*pixels/25 - [0 scaleBarSize*scale],0.75*pixels/25 + [0 0],...
    'Color','w','LineWidth',3);
scaleText = text(24.25*pixels(1)/25 - (scaleBarSize/2)*scale,pixels(2)/25,[num2str(scaleBarSize),' mm'],...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Color','w','FontName','Arial','FontSize',24);

% Plot time (s)
timeText = text(24.25*pixels(1)/25,24.25*pixels(2)/25,datestr(seconds(0),'HH:MM:SS'),...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'Color','w','FontName','Arial','FontSize',24);
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],...
    'XTick',[],'YTick',[],'Units','inches')
set(gcf,'Toolbar','none')

% Plot worm
p = plot(wormX(1,:),wormY(1,:),'k','LineWidth',2);
s = scatter(wormX(1,1),wormY(1,1),150,'k.');

for i = 1:length(wormX)
    p.XData = wormX(i,:); p.YData = wormY(i,:);
    s.XData = wormX(i,1); s.YData = wormY(i,1);

    if any(all(wormTrack(i,1) >= encounter.enter(indEncounter) & wormTrack(i,1) <= encounter.exit(indEncounter),2))
        p.Color = [0 0 1]; s.CData = [0 0 1];
    else
        p.Color = [0 0 0]; s.CData = [0 0 0];
    end

    set(timeText,'String',datestr(seconds(wormTrack(i,2)),'HH:MM:SS'));

    % Save frame
    frame = getframe(gca);
    writeVideo(v,frame);
end

% Close video
close(v)

%% Video S5 - Sensing v. Non-Sensing 3D plot

% code in figS11.m

%% Video S6 - Tracking video w/ labeled encounters

wormNum = 360;
load(fullfile(path,'foragingConcentration','experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
indInfo = find(cellfun(@(w) any(w == wormNum),info.wormNum));
indEncounter = find(strcmp(encounter.expName,'foragingConcentration') & encounter.wormNum == wormNum);
indData = data.wormNum == wormNum;
video = getVideoInfo(info.videoFileName{indInfo});
videoLength = 120; % seconds
numSteps = 30*videoLength;

wormData = data(indData,:);
wormX = round(wormData.xPosition);
wormY = round(wormData.yPosition);
wormFrame = wormData.frameNum;
wormTime = wormData.timeOffset;
edges = round(linspace(1,max(wormFrame),numSteps));

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
encounterColors = encounter.exploitPosterior(indEncounter).*exploitColor + ...
            (1-encounter.exploitPosterior(indEncounter)).*encounter.sensePosterior(indEncounter).*sampleColor + ...
            (1-encounter.exploitPosterior(indEncounter)).*(1-encounter.sensePosterior(indEncounter)).*searchOnColor;

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
lawnOD600 = info.lawnClosestOD600{indInfo}.*info.lawnMask{indInfo};
lawnBorder = ones(size(lawnOD600));
for i = 1:length(uniqueOD600)
    lawnBorder(lawnOD600 == uniqueOD600(i)) = colorValue(i);
end
lawnImage = uint8(255.*lawnBorder);
arenaMask = uint8(info.arenaMask{indInfo});

v = VideoWriter([saveDir,'VideoS6'],'MPEG-4');
open(v);

figure('Position',[0 0 560 560])
displayImage = flip(arenaMask.*lawnImage);
displayImage(displayImage == 0) = 34; % Keynote gray
displayImage = repmat(displayImage,1,1,3);
h = imshow(displayImage);
hold on
scaleBarSize = 5; % mm
scale = info.scale(indInfo);
pixels = info.pixels{indInfo};

% Plot scale bar
scaleBar = plot(24.25*pixels/25 - [0 scaleBarSize*scale],0.75*pixels/25 + [0 0],...
    'Color','w','LineWidth',3);
scaleText = text(24.25*pixels(1)/25 - (scaleBarSize/2)*scale,pixels(2)/25,[num2str(scaleBarSize),' mm'],...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Color','w','FontName','Arial','FontSize',24);

% Plot time (s)
timeText = text(24.25*pixels(1)/25,24.25*pixels(2)/25,datestr(seconds(0),'HH:MM:SS'),...
    'HorizontalAlignment','right','VerticalAlignment','top',...
    'Color','w','FontName','Arial','FontSize',24);
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],...
    'XTick',[],'YTick',[],'Units','inches')
set(gcf,'Toolbar','none')

% Plot state name
stateText = text(0.75*pixels(1)/25,24.25*pixels(2)/25,'search',...
    'HorizontalAlignment','left','VerticalAlignment','top',...
    'Color','w','FontName','Arial','FontSize',36,'FontWeight','bold');

% Set track filter
filterRadius = 2;
structure = strel('disk',filterRadius); structure = structure.Neighborhood;
[filterX,filterY] = ind2sub(size(structure),find(structure));
filterX = filterX - filterRadius; filterY = filterY - filterRadius;

% Plot worm
for i = 1:length(wormX)
    ind = find(all(wormFrame(i) >= encounter.enter(indEncounter) & ...
        wormFrame(i) <= encounter.exit(indEncounter),2));
    if ~any(isnan([wormX(i),wormY(i)]))
        position = sub2ind(pixels,wormY(i)+filterY,wormX(i)+filterX) + [0 1 2].*prod(pixels);
        if ~isempty(ind)
            displayImage(position) = repmat(255*encounterColors(ind,:),length(filterY),1);
            if encounter.sensePosterior(indEncounter(ind)) < 0.5
                set(stateText,'String','search','Color',encounterColors(ind,:));
            elseif encounter.exploitPosterior(indEncounter(ind)) >= 0.5
                set(stateText,'String','exploit','Color',encounterColors(ind,:));
            else
                set(stateText,'String','sample','Color',encounterColors(ind,:));
            end
        else
            displayImage(position) = zeros(length(filterY),3);
            set(stateText,'String','search','Color','w');
        end
        set(h,'CData',displayImage);
    end
    set(timeText,'String',datestr(seconds(wormTime(i)),'HH:MM:SS'));

    % Save frame
    if ismember(i,edges)
        frame = getframe(gca);
        writeVideo(v,frame);
    end
end

% Close video
close(v)