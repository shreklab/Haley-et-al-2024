%% Load Data

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
saveDir = [path,'figures\FigS5\'];
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,'trajectory.mat'),'trajectory');
% data = readtable([path,expName,'\videos\22-03-18\data\midpoint\2022-03-18_13-04-43_2.csv']);

%% Get conditions to analyze

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));

conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 0.5 & ...
    ~(strcmp(GID.OD600Label,'0.10') & GID.growthCondition == 48));
numGroups = length(conditionG);

%% Get ids for worms used in this model

wormNums = unique(encounter.wormNum(ismember(G,conditionG) & ~encounter.exclude));

% Remove worms tracked for less than 75% of the video
framesTracked = arrayfun(@(w) sum(~data.noTrack(data.wormNum == w))/...
   sum(info.numFrames(info.plateNum == unique(data.plateNum(data.wormNum == w)))), wormNums);
wormNums = wormNums(framesTracked >= 0.75);
numWorms = length(wormNums);

% Get condition id of each worm
wormGroup = arrayfun(@(w) unique(G(encounter.wormNum == w & ismember(G,conditionG))),wormNums);
[~,~,wormGroupInd] = unique(wormGroup); % continuous group idviolin

% Get indices of encounters labeled exploit, sample, or searchOn
indEncounter = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude;

%% FigS5A - Quantification of Deceleration (Average)

winTime = [60 240]; % -60 to +180 seconds from enter/exit event
timeWindow = -winTime(1):0.25:winTime(2);
indTime = find(timeWindow >= -20 & timeWindow <= 40);
timeWindow = timeWindow(indTime);

indTrajectory = join(trajectory.enterInfo,encounter(strcmp(encounter.expName,expName),:),'Keys',...
    {'expNum','plateNum','videoNum','wormNum','id'});
indTrajectory = ismember(indTrajectory.wormNum,wormNums) & ~indTrajectory.exclude;

enterVelocity = trajectory.enterVelocity(indTrajectory,indTime);
enterVelocityMean = mean(enterVelocity,'omitnan');
indDeceleration = find(timeWindow == -1.5 | timeWindow == 6.5);

figure; hold on;
patch([0 0 40 40],[0 400 400 0], 0.4055.*[1 1 1],'EdgeAlpha',0,'FaceAlpha',0.5)
plot(timeWindow,enterVelocityMean,'k','LineWidth',1)
scatter(timeWindow(indDeceleration(1)),enterVelocityMean(indDeceleration(1)),'g')
scatter(timeWindow(indDeceleration(2)),enterVelocityMean(indDeceleration(2)),'r')
xline(-1.5,'g'); xline(6.5,'r')
xlim([-20 40]); xticks(-20:20:40); xticklabels(-20:20:40); ylim([150 250])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
xlabel('time to patch (s)','FontSize',8); ylabel('velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 1.5])
exportgraphics(gcf,[saveDir,'figS5a.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S5B - Quantification of Deceleration (Example 1)

wormNum = 360;
encounterNum = 6;
indInfo = find(cellfun(@(w) any(w == wormNum),info.wormNum));
indEncounter = find(strcmp(encounter.expName,expName) & encounter.wormNum == wormNum & encounter.id > 0);
indData = data.wormNum == wormNum;
time = data.time(indData);
velocitySmooth = data.velocitySmooth(indData);

% Compute deceleration at slowdown via linear fit from -1.5 to +6.5
% seconds from patch entry
indDecelerate = time >= encounter.timeEnter(indEncounter(encounterNum)) - 1.5 & ...
    time <= encounter.timeEnter(indEncounter(encounterNum)) + 6.5;
dT = time(indDecelerate) - mean(time(indDecelerate),'omitnan');
dV = velocitySmooth(indDecelerate) - mean(velocitySmooth(indDecelerate),'omitnan');
decelerate = sum(dT.*dV,'omitnan')./sum(dT.^2,'omitnan');

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

figure; hold on;
patch(encounter.timeEnter(indEncounter(encounterNum)).*[1 1 1.5 1.5],[0 400 400 0],...
    colorValue.*[1 1 1],'EdgeAlpha',0,'FaceAlpha',0.5)
plot(time,velocitySmooth,'k','LineWidth',1)
xline(time(find(indDecelerate,1,'first')),'g')
xline(time(find(indDecelerate,1,'last')),'r')
f = fit(time(indDecelerate),velocitySmooth(indDecelerate),'poly1');
plot(time(indDecelerate),f(time(indDecelerate)),'b','LineWidth',1)
xlim(encounter.timeEnter(indEncounter(encounterNum)) + [-20 40])
xticks(encounter.timeEnter(indEncounter(encounterNum)) + [-20:20:40]);
xticklabels(-20:20:40); ylim([0 400])
text(mean(time(indDecelerate),'omitnan'),mean(velocitySmooth(indDecelerate),'omitnan'),...
    [num2str(decelerate),' (μm/s^2)'],'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Color','b','FontName','Arial','FontSize',8);
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
xlabel('time to patch (s)','FontSize',8); ylabel('velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 1.5])
exportgraphics(gcf,[saveDir,'figS5b.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S5B - Quantification of Deceleration (Example 2)

encounterNum = 20;
indInfo = find(cellfun(@(w) any(w == wormNum),info.wormNum));
indEncounter = find(strcmp(encounter.expName,expName) & encounter.wormNum == wormNum & encounter.id > 0);
indData = data.wormNum == wormNum;
time = data.time(indData);
velocitySmooth = data.velocitySmooth(indData);

% Compute deceleration at slowdown via linear fit from -1.5 to +6.5
% seconds from patch entry
indDecelerate = time >= encounter.timeEnter(indEncounter(encounterNum)) - 1.5 & ...
    time <= encounter.timeEnter(indEncounter(encounterNum)) + 6.5;
dT = time(indDecelerate) - mean(time(indDecelerate),'omitnan');
dV = velocitySmooth(indDecelerate) - mean(velocitySmooth(indDecelerate),'omitnan');
decelerate = sum(dT.*dV,'omitnan')./sum(dT.^2,'omitnan');

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

figure; hold on;
patch(encounter.timeEnter(indEncounter(encounterNum)).*[1 1 1.5 1.5],[0 400 400 0],...
    colorValue.*[1 1 1],'EdgeAlpha',0,'FaceAlpha',0.5)
plot(time,velocitySmooth,'k','LineWidth',1)
xline(time(find(indDecelerate,1,'first')),'g')
xline(time(find(indDecelerate,1,'last')),'r')
f = fit(time(indDecelerate),velocitySmooth(indDecelerate),'poly1');
plot(time(indDecelerate),f(time(indDecelerate)),'b','LineWidth',1)
xlim(encounter.timeEnter(indEncounter(encounterNum)) + [-20 40])
xticks(encounter.timeEnter(indEncounter(encounterNum)) + [-20:20:40]);
xticklabels(-20:20:40); ylim([0 400])
text(mean(time(indDecelerate),'omitnan'),mean(velocitySmooth(indDecelerate),'omitnan'),...
    [num2str(decelerate),' (μm/s^2)'],'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Color','b','FontName','Arial','FontSize',8);
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
xlabel('time to patch (s)','FontSize',8); ylabel('velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 1.5])
exportgraphics(gcf,[saveDir,'figS5c.pdf'],'BackgroundColor','none','ContentType','vector')