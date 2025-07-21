function [] = plotTracks(info,data,wormNum,metric,saveDir,bodyPart,overwrite)
%[] = PLOTTRACKS(info,data,wormNum,metric,saveDir,bodyPart,overwrite)
%
%   PLOTTRACKS writes an image to file showing the tracks of an animal
%   overlaid onto an image of the lawns and arena. Tracks are colored based
%   on a given metric (e.g., time, velocity, path angle).
%
%   INPUTS:
%       - info [struct]: structure containing the metadata for all worm(s);
%           created with GETFORAGINGINFO
%       - data [struct]: structure containing the location of worm(s) at 
%           every timepoint; created with ANALYZEFORAGING
%       - wormNum [double]: identifying worm # to plot
%       - metric [str]: name of metric as found in the data structure
%           (e.g., 'timeOffset', 'velocitySmooth', or 'pathAngle').
%       - saveDir [str]: path name to save the image file
%       - bodyPart [str]: part of the animal to plot in the event that
%           multiple body segments were tracked (e.g.,'head','midpoint', or
%           'tail')
%       - overwrite [str]: overwrite the file (i.e., 'Yes' or 'No'); if
%           'No', this function does not plot
%
%   Written 5/8/2023 by Jess Haley in MATLAB R2023a.
%
%   See also ANALYZEFORAGING, GETFORAGINGINFO, PLOTTRACKS.

warning('off','MATLAB:print:ContentTypeImageSuggested')
maxTime = 1; % hours

% Get plate and video number
plateNum = data.plateNum(find(data.wormNum == wormNum,1));
videoNum = data.videoNum(find(data.wormNum == wormNum,1));
indInfo = info.plateNum == plateNum & info.videoNum == videoNum;

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

% Check if folder exists to save figures
figureFolder = [saveDir,datestr(info.timeRecord(indInfo),'yy-mm-dd'),...
    filesep,'plot',filesep,bodyPart,filesep,metric,filesep];
if ~exist(figureFolder,'dir')
    mkdir(figureFolder)
end

% Check if file exists and return out of function if not overwriting
figureFile = [figureFolder,info.videoFileName{indInfo}(1:end-4),...
    '_',num2str(wormNum,'%.4d'),'.pdf'];
if exist(figureFile,'file') && strcmp(overwrite,'No')
    return
end

% Get Data
% colorData = data.(metric)(data.wormNum == wormNum);
% xData = data.xPosition(data.wormNum == wormNum);
% yData = data.yPosition(data.wormNum == wormNum);
% turnData = data.turn(data.wormNum == wormNum);
colorData = data.(metric)(data.wormNum == wormNum & data.timeOffset <= 3600*maxTime);
xData = data.xPosition(data.wormNum == wormNum & data.timeOffset <= 3600*maxTime);
yData = data.yPosition(data.wormNum == wormNum & data.timeOffset <= 3600*maxTime);
turnData = data.turn(data.wormNum == wormNum & data.timeOffset <= 3600*maxTime);
% lawnGrowth = info.lawnGrowth(indInfo);
lawnOD600 = info.lawnClosestOD600{indInfo}.*info.lawnMask{indInfo};
lawnBorder = ones(size(lawnOD600));
for i = 1:length(uniqueOD600)
    % lawnBorder(lawnOD600 == uniqueOD600(i)) = 1 - (log(borderValue(i))*(0.6/8) + 0.35);
    lawnBorder(lawnOD600 == uniqueOD600(i)) = colorValue(i);
end
lawnImage = uint8(255.*lawnBorder);
% lawnImage = uint8(sqrt(2*days(lawnGrowth).*lawnOD600).*-40+255);
arenaMask = uint8(info.arenaMask{indInfo});

% Set color bins
if strcmp(metric,'time') || strcmp(metric,'pathAngle') || strcmp(metric,'timeOffset')
    colorData = abs(colorData);
    bins = ceil(linspace(1,max(colorData),60));
elseif contains(metric,'velocity') > 0
    colorData(colorData > 500) = NaN;
    bins = round(linspace(0,400,59));
    bins(end+1) = 500;
end
c = jet(length(bins));

% Set figure properties
figure('Position',[0 0 560 560])
displayImage = flip(arenaMask.*lawnImage);
displayImage(displayImage == 0) = 34; % Keynote gray
imshow(displayImage)
hold on
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],...
    'XTick',[],'YTick',[])

% Plot worm position bin-by-bin
if strcmp(metric,'time') || strcmp(metric,'timeOffset')
    for j = 1:length(bins)-1
        ind = find(colorData >= bins(j));
        plot(xData(ind),yData(ind),'Color',c(j,:),'LineWidth',1);
    end
elseif contains(metric,'velocity') || strcmp(metric,'pathAngle')
    for j = 2:length(xData)
        if isnan(colorData(j))
            plot(xData(j+[-1,0]),yData(j+[-1,0]),'Color','k','LineWidth',1)
        else
            thisColor = c(find(bins >= colorData(j),1),:);
            plot(xData(j+[-1,0]),yData(j+[-1,0]),...
                'Color',thisColor,'LineWidth',1)
        end
    end
end
if contains(metric,'pathAngle')
    scatter(xData(turnData),yData(turnData),'k')
end

% Export figure
set(gcf,'Toolbar','none')
exportgraphics(gca,figureFile,'BackgroundColor','none','ContentType','vector')
% figureFile = [figureFolder,info.videoFileName{indInfo}(1:end-4),...
%     '_',num2str(wormNum,'%.4d'),'.png'];
% exportgraphics(gca,figureFile)
close(gcf);

end