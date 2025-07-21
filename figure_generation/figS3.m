%% Load Data

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
saveDir = [path,'figures\FigS3\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
permuteColor = [204 121 167]./255; % pink

%% Get condition ids

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume','growthCondition',...
    'lawnOD600','peptone'}));
conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 0.5 & ...
    GID.lawnOD600 == 10);
numGroups = length(conditionG);

%% Get average estimated amplitude of each condition and normalize to OD = 10 (0.5 uL)

borderAmp = splitapply(@(X) mean(X,'omitnan'),encounter.borderAmplitude,G);
relativeBorder = borderAmp(conditionG);
borderAmp10 = borderAmp(strcmp(GID.expName,'foragingConcentration') & ...
    GID.lawnOD600 == 10.00 & GID.lawnVolume == 0.5);
borderAmp0 = 1e-2; % assign 0 to 0.01
relativeBorder = 10.*relativeBorder./borderAmp10;
relativeBorder(GID.lawnOD600(conditionG) == 1e-10) = borderAmp0;

% Assign each border value to a color
colorValue = min(1 - (log(relativeBorder)*0.09 + 0.4),0.8);

%% Get ids for worms used in this figure

wormNums = unique(encounter.wormNum(ismember(G,conditionG) & ~encounter.exclude));

% Remove worms tracked for less than 75% of the video
framesTracked = arrayfun(@(w) sum(~data.noTrack(data.wormNum == w))/...
   sum(info.numFrames(info.plateNum == unique(data.plateNum(data.wormNum == w)))), wormNums);
wormNums = wormNums(framesTracked >= 0.75);
numWorms = length(wormNums);

% Get condition id of each worm
wormGroup = arrayfun(@(w) unique(G(encounter.wormNum == w & ismember(G,conditionG))),wormNums);
[~,~,wormGroupInd] = unique(wormGroup); % continuous group id

% Get indices of encounters labeled exploit, sample, or searchOn
indEncounter = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn'));

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%% Figure S3A - Encounter Duration Gaussian Mixture Model

patchDuration = log10(encounter.duration(indEncounter)./60);
mdl = fitgmdist(patchDuration,2,'Replicates',10000);
durationPosterior = posterior(mdl,patchDuration);
[~,longCluster] = max(mdl.mu);
[~,shortCluster] = min(mdl.mu);
patchLong = durationPosterior(:,longCluster);
patchShort = durationPosterior(:,shortCluster);
x = logspace(log10(1/30),log10(60),1000);
yLong = normpdf(log10(x),mdl.mu(longCluster),sqrt(mdl.Sigma(longCluster)));
yShort = normpdf(log10(x),mdl.mu(shortCluster),sqrt(mdl.Sigma(shortCluster)));
yModel = pdf(mdl,log10(x)');

[10.^mdl.mu(shortCluster),10.^mdl.Sigma(shortCluster)]
[10.^mdl.mu(longCluster),10.^mdl.Sigma(longCluster)]

figure; hold on;
histogram(patchDuration,-1.5:0.1:1.8,...
    'Normalization','pdf','FaceColor',colorValue.*[1 1 1])
%plot(log10(x),yModel,'k','LineWidth',2)
plot(log10(x),yShort.*mdl.ComponentProportion(shortCluster),'Color',sampleColor,'LineWidth',2)
plot(log10(x),yLong.*mdl.ComponentProportion(longCluster),'Color',exploitColor,'LineWidth',2)
xticks(log10([2 10 60 600 3600]./60))
xticklabels({'1/30','1/6','1','10','60'});
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','YTick',[]); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('encounter duration (min)','FontSize',8); ylabel('pdf','FontSize',8)
exportgraphics(gca,[saveDir,'figS3a.pdf'])

%% Figure S3B - Posterior Probabilities of GMM Clusters

figure; hold on;
plot(sort(patchShort,'descend'),'LineWidth',2,'Color',sampleColor)
plot(sort(patchLong,'ascend'),'LineWidth',2,'Color',exploitColor)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xticks([0 length(patchDuration)]); yticks([0 1])
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('# encounters','FontSize',8); ylabel('p(cluster)','FontSize',8)
exportgraphics(gca,[saveDir,'figS3b.pdf'])

%% Figure S3C - Encounters labeled for Censorship

rng(1);
wormNums = wormNums(randperm(numWorms));
wormNums = [360;wormNums(wormNums ~= 360)];
numSec = 3601;
colorMap = uint8(255.*[permuteColor;colorValue.*ones(1,3)]);
durationMap = linspace(-1.5,1.8,100);

% Get patch colors
patchDurationColor = uint8(255.*ones(numWorms,numSec,3));
for i = 1:numWorms
    ind = find(encounter.wormNum == wormNums(i) & indEncounter);
    for j = 1:length(ind)
        timeWindow = round(encounter.timeEnter(ind(j))):round(encounter.timeExit(ind(j)));
        censor = encounter.censorEnter(ind(j)) | encounter.censorExit(ind(j));
        [~,indColor] = min(abs(log10(encounter.duration(ind(j))./60) - durationMap));
        patchDurationColor(i,timeWindow+1,:) = repmat(colorMap(censor+1,:),length(timeWindow),1);
    end
end
patchDurationColor = imresize(patchDurationColor,[numWorms*round((0.75/2.5)*numSec/numWorms) numSec],...
    'method','nearest','Colormap','original');

figure;
imshow(patchDurationColor)
ax = get(gca);
set(gca,'Units','inches','Visible','on','YTick',ax.YLim,'YTickLabel',[1 numWorms],...
    'XTick',0.5:600:3600.5,'XTickLabel',0:10:60,'FontName','Arial','FontSize',8,...
    'GridLineWidth',0.75)
ylabel('worm #','FontSize',8); xlabel('time (min)','FontSize',8)
set(gca,'Position',[1 1 2.5 0.75])
exportgraphics(gca,[saveDir,'figS3c.pdf'])

%% Figure S3D - Encounter Duration (Uncensored Only)

indUncensored =  ~encounter.censorEnter & ~encounter.censorExit;
numCensorEnter = sum(indEncounter & encounter.censorEnter)/numWorms
numCensorExit = sum(indEncounter & encounter.censorExit)/numWorms
pCensor = sum(indEncounter & ~indUncensored)/sum(indEncounter)

patchDurationUncensor = log10(encounter.duration(indEncounter & indUncensored)./60);
mdlUncensor = fitgmdist(patchDurationUncensor,2,'Replicates',10000);
durationPosteriorUncensor = posterior(mdlUncensor,patchDurationUncensor);
[~,longCluster] = max(mdlUncensor.mu);
[~,shortCluster] = min(mdlUncensor.mu);
patchLong = durationPosteriorUncensor(:,longCluster);
patchShort = durationPosteriorUncensor(:,shortCluster);
x = logspace(log10(1/30),log10(60),1000);
yLong = normpdf(log10(x),mdlUncensor.mu(longCluster),sqrt(mdlUncensor.Sigma(longCluster)));
yShort = normpdf(log10(x),mdlUncensor.mu(shortCluster),sqrt(mdlUncensor.Sigma(shortCluster)));
yModel = pdf(mdlUncensor,log10(x)');

[10.^mdlUncensor.mu(shortCluster),10.^mdlUncensor.Sigma(shortCluster)]
[10.^mdlUncensor.mu(longCluster),10.^mdlUncensor.Sigma(longCluster)]

figure; hold on;
histogram(patchDurationUncensor,-1.5:0.1:1.8,...
    'Normalization','pdf','FaceColor',permuteColor)
%plot(log10(x),yModel,'k','LineWidth',2)
plot(log10(x),yShort.*mdlUncensor.ComponentProportion(shortCluster),'Color',sampleColor,'LineWidth',2)
plot(log10(x),yLong.*mdlUncensor.ComponentProportion(longCluster),'Color',exploitColor,'LineWidth',2)
xticks(log10([2 10 60 600 3600]./60))
xticklabels({'1/30','1/6','1','10','60'});
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','YTick',[]); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('encounter duration (min)','FontSize',8); ylabel('pdf','FontSize',8)
exportgraphics(gca,[saveDir,'figS3d.pdf'])

%% Figure S3E - Posterior Probabilities of GMM Clusters (Uncensored Only)

figure; hold on;
plot(sort(patchShort,'descend'),'LineWidth',2,'Color',sampleColor)
plot(sort(patchLong,'ascend'),'LineWidth',2,'Color',exploitColor)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xticks([0 length(patchDurationUncensor)]); yticks([0 1]); xlim([0 length(patchDurationUncensor)])
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('# encounters','FontSize',8); ylabel('p(cluster)','FontSize',8)
exportgraphics(gca,[saveDir,'figS3e.pdf'])

%% Figure S3F - Time Enter (Uncensored Only)

timeEnter = encounter.timeEnter(indEncounter & indUncensored)./60;
pTimeEnter = ranksum(timeEnter(patchLong < 0.5),timeEnter(patchLong >= 0.5))
shortEnter = median(timeEnter(patchLong < 0.5))
longEnter = median(timeEnter(patchLong >= 0.5))

figure; hold on
v = Violin({timeEnter(patchLong < 0.5)},1,'ViolinColor',{sampleColor},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
v = Violin({timeEnter(patchLong >= 0.5)},2,'ViolinColor',{exploitColor},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
set(gca,'Units','inches','XTick',1:2,'XTickLabel',{'short','long'},...
    'FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975*35/65 0.75])
text(1.5,59,'***','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',8,'FontName','Arial')
plot([1 2],[65 65],'k','LineWidth',1)
xlim([0.4 2.6]); ylim([0 75]); yticks(0:20:60)
ylabel('encounter start (min)','FontSize',8)
exportgraphics(gca,[saveDir,'figS3f.pdf'])

%% Figure S3G - Comparison between Models

durationPosterior = posterior(mdl,patchDuration);
[~,longCluster] = max(mdlUncensor.mu);
patchLong = durationPosterior(:,longCluster);

durationPosteriorUncensor = posterior(mdlUncensor,patchDuration);
[~,longCluster] = max(mdl.mu);
patchLongUncensor = durationPosteriorUncensor(:,longCluster);
clusterChange = sum(all([patchLong,patchLongUncensor] > 0.5,2)-any([patchLong,patchLongUncensor] > 0.5,2)~=0)

figure; hold on;
plot([0 1],[0,1],'k','LineWidth',1)
patch([0 0.5 0.5 0],[0 0 0.5 0.5],'k','FaceAlpha',0.3,'EdgeAlpha',0)
patch([0.5 1 1 0.5],[0.5 0.5 1 1],'k','FaceAlpha',0.3,'EdgeAlpha',0)
plot([patchLong,patchLong]',[patchLong,patchLongUncensor]','Color',permuteColor,'LineWidth',1)
scatter(patchLong,patchLong,'filled','MarkerFaceColor','k');
scatter(patchLong,patchLongUncensor,'filled','MarkerFaceColor',permuteColor);
ylim([0 1])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.75 0.75])
xlabel('p(long)_{all}','FontSize',8); ylabel('p(long)','FontSize',8)
exportgraphics(gca,[saveDir,'figS3g.pdf'])