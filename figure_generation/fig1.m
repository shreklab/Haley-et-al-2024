%% Load data into workspace

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
load(fullfile(path,expName,'trajectory.mat'),'trajectory');
load(fullfile(path,expName,'distance.mat'),'distance');
load(fullfile(path,expName,'permutePatchLocation.mat'),'onPatchShuffled');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
thresh = readtable([path,'foragingMini\encounterThresholds.csv']);
saveDir = [path,'figures\Fig1\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
searchOffColor = [240 228 66]./255; % yellow
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
    ~encounter.exclude & ~strcmp(encounter.label,'') & ~strcmp(encounter.label,'searchOff');

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%% Figure 1A - Example Worm Trace (Time)

wormNum = 360;
% plotTracks(info,data,wormNum,'timeOffset',saveDir,bodyPart,'yes');

%% Figure 1B - Example Worm Distance to Patch Edge

wormNum = 360;
indWorm = find(data.wormNum == wormNum);
indInfoWorm = find(info.plateNum == data.plateNum(indWorm(1)));
indEncounterWorm = find(encounter.wormNum == wormNum & strcmp(encounter.expName,expName) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn')));
scale = info.scale(indInfoWorm);
arenaEdge = bwdist(info.lawnMask{indInfoWorm}).*bwperim(info.arenaMask{indInfoWorm})./scale;
arenaEdge = -round(max(arenaEdge(arenaEdge>0)),1);
patchMiddle = round(max(bwdist(~info.lawnMask{indInfoWorm}),[],'all')./scale,1);
distanceLawnEdge = data.distanceLawnEdge(indWorm)./scale;
distanceLawnEdge(data.noTrack(indWorm) == 1) = NaN;

figure('Position',[360 278 840 420]);
hold on
yline(0,'k--','LineWidth',0.75);
for i = 1:length(indEncounterWorm)
    patch([encounter.timeEnter(indEncounterWorm(i)).*[1 1] ...
        encounter.timeExit(indEncounterWorm(i)).*[1 1]],...
        [arenaEdge patchMiddle patchMiddle arenaEdge],colorValue.*[1 1 1],...
        'EdgeColor','none')
end
plot(data.timeOffset(indWorm),distanceLawnEdge,'k','LineWidth',0.5)
xlim([0 3600]); xticks(0:600:3600); xticklabels(0:10:60); 
ylim([arenaEdge patchMiddle]); yticks(ceil(arenaEdge):1); 
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2.5 0.75])
xlabel('time (min)','FontSize',8);ylabel('dist. to patch (mm)','FontSize',8)
exportgraphics(gca,[saveDir,'fig1b.pdf'])

%% Figure 1C - Encounters for All OD600 10 Worms

rng(1);
wormNums = wormNums(randperm(numWorms));
wormNums = [360;wormNums(wormNums ~= 360)];
numSec = 3601;
colorMap = uint8(cell2mat(arrayfun(@(c1,c2) linspace(c1,c2,100)',...
    255.*searchOffColor,255.*exploitColor,'UniformOutput',false)));
% colorMap = flipud(uint8(cell2mat(arrayfun(@(c1,c2) log10(logspace(c2,c1,100))',...
%     255.*searchOffColor,max(255.*exploitColor,1),'UniformOutput',false))));
durationMap = linspace(-1.5,1.8,100);
%durationMap = logspace(log10(1/30),log10(60),100);
%durationMap = linspace(1/30,60,100);

% Get patch colors
patchColors = uint8(255.*ones(numWorms,numSec));
onPatchData = false(numWorms,numSec);
patchDurationColor = uint8(255.*ones(numWorms,numSec,3));
for i = 1:numWorms
    ind = find(encounter.wormNum == wormNums(i) & indEncounter);
    onPatch = zeros(1,numSec);
    onPatch(floor(encounter.timeEnter(ind))+1) = 1;
    onPatch(floor(encounter.timeExit(ind))+2) = -1;
    patchColors(i,cumsum(onPatch) > 0) = 255*colorValue;
    for j = 1:length(ind)
        timeWindow = round(encounter.timeEnter(ind(j))):round(encounter.timeExit(ind(j)));
        [~,indColor] = min(abs(log10(encounter.duration(ind(j))./60) - durationMap));
        %[~,indColor] = min(abs(encounter.duration(ind(j))/60 - durationMap));
        patchDurationColor(i,timeWindow+1,:) = repmat(colorMap(indColor,:),length(timeWindow),1);
        onPatchData(i,timeWindow+1) = 1;
    end
end
patchColors = imresize(patchColors,[numWorms*round((0.75/2.5)*numSec/numWorms) numSec],...
    'method','nearest','Colormap','original');
patchDurationColor = imresize(patchDurationColor,[numWorms*round((0.75/2.5)*numSec/numWorms) numSec],...
    'method','nearest','Colormap','original');

figure;
imshow(patchDurationColor)
ax = get(gca);
set(gca,'Units','inches','Visible','on','YTick',ax.YLim,'YTickLabel',[1 numWorms],...
    'XTick',0.5:600:3600.5,'XTickLabel',0:10:60,'FontName','Arial','FontSize',8,...
    'GridLineWidth',0.75)
colormap(colorMap); c = colorbar(gca); clim(prctile(durationMap,[0 100]));
c.Ticks = log10([2 10 60 600 3600]./60); c.TickLabels = {'1/30','1/6','1','10','60'};
c.Label.String = 'duration (min)'; c.Label.Rotation = -90;
ylabel('worm #','FontSize',8); xlabel('time (min)','FontSize',8)
set(gca,'Position',[1 1 2.5 0.75])
exportgraphics(gca,[saveDir,'fig1c.pdf'])

%% Figure 1D - # Patch Visits (total + unique)

figure;
numEncounters = arrayfun(@(w) sum(indEncounter & encounter.wormNum == w), wormNums);
uniqueEncounters = arrayfun(@(w) length(unique(encounter.lawnID(indEncounter & encounter.wormNum == w))),wormNums);
subplot(4,4,[1,5]); hold on
v = Violin({numEncounters},1,'ViolinColor',{colorValue.*[1 1 1]},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
v = Violin({uniqueEncounters},2,'ViolinColor',{colorValue.*[1 1 1]},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
set(gca,'Units','inches','XTick',1:2,'XTickLabel',{'total','unique'},...
    'FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975*35/65 0.75])
xlim([0.4 2.6]); ylim([0 30]); yticks(0:10:30)
ylabel('# encounters','FontSize',8)
exportgraphics(gca,[saveDir,'fig1d.pdf'])

% Source Data
summaryTableD = table([numWorms;numWorms],[mean(numEncounters);mean(uniqueEncounters)],...
    [std(numEncounters);std(uniqueEncounters)],[median(numEncounters);median(uniqueEncounters)],...
    [prctile(numEncounters,2.5);prctile(uniqueEncounters,2.5)],[prctile(numEncounters,97.5);prctile(uniqueEncounters,97.5)],...
    'VariableNames',{'# worms','mean','std','median','2.5%','97.5%'},'RowNames',{'total','unique'},...
    'DimensionNames',{'encounter type','Variables'})
sourceTableD = table((1:numWorms)',numEncounters,uniqueEncounters,'VariableNames',...
    {'worm #','total encounters','unique encounters'});
writeSourceData([saveDir,'fig1_data1.xlsx'],{'Figure 1 - Source Data 1';'Figure 1D - Number of total and unique patch encounters'},...
    {'Summary Statistics';'Source Data'},{summaryTableD;sourceTableD},{[];[]})

%% Figure 1E,F - Encounter Duration Histogram (OD600 = 10) & Time Encounter

indUncensored =  ~encounter.censorEnter & ~encounter.censorExit;
numCensorEnter = sum(indEncounter & encounter.censorEnter)/numWorms
numCensorExit = sum(indEncounter & encounter.censorExit)/numWorms
pCensor = sum(indEncounter & ~indUncensored)/sum(indEncounter)

% patchDuration = log10(encounter.duration(indEncounter & indUncensored)./60);
patchDuration = log10(encounter.duration(indEncounter)./60);
mdl = fitgmdist(patchDuration,2);
patchLong = posterior(mdl,patchDuration);
probVar = mean(prod(patchLong,2))
[~,exploitCluster] = max(mdl.mu);
patchLong = patchLong(:,exploitCluster);

subplot(4,4,[2:4,6:8]); hold on;
h1 = histogram(patchDuration(patchLong >= 0.5),-1.5:0.1:1.8,...
    'Normalization','count','FaceColor',exploitColor);
h2 = histogram(patchDuration(patchLong < 0.5),-1.5:0.1:1.8,...
    'Normalization','count','FaceColor',sampleColor);
xticks(log10([2 10 60 600 3600]./60))
xticklabels({'1/30','1/6','1','10','60'});
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('encounter duration (min)','FontSize',8); ylabel('# encounters','FontSize',8)
exportgraphics(gca,[saveDir,'fig1e.pdf'])

% Source Data
summaryTableE = table([numWorms;numWorms],[sum(patchLong < 0.5);sum(patchLong >= 0.5)],...
    mdl.mu,reshape(mdl.Sigma,2,1),reshape(mdl.ComponentProportion,2,1),...
    'VariableNames',{'# worms','# encounters','mu','sigma','proportion'},'RowNames',{'short','long'},...
    'DimensionNames',{'encounter type','Variables'});
summaryNoteE = {'* mu, sigma, and proportion define the mean, covariance, and component proportion of each gaussian';...
    '* mu and sigma are in units of log10(encounter duration in minutes)'};
histogramTableE = table((1:h1.NumBins)',10.^(h1.BinEdges(1:end-1) + diff(h1.BinEdges(1:2))/2)',...
    (h1.BinEdges(1:end-1) + diff(h1.BinEdges(1:2))/2)',h1.Values' + h2.Values',...
    [repmat({'short'},find(h2.Values > 0,1,'last'),1);repmat({'long'},h1.NumBins - find(h2.Values > 0,1,'last'),1)],...
    'VariableNames',{'bin #','duration (min)','log10(duration)','# encounters','encounter type'})
[~,wormList] = ismember(encounter.wormNum(indEncounter),wormNums);
encounterLabel = {'short','long'}; encounterLabel = encounterLabel(1+(patchLong >= 0.5))';
encounterList = arrayfun(@(w) 1:sum(encounter.wormNum(indEncounter)==w),unique(encounter.wormNum(indEncounter)),'UniformOutput',false);
encounterList = [encounterList{:}]';
sourceTableE = table(wormList,encounterList,encounter.timeEnter(indEncounter)./60,...
    10.^patchDuration,patchDuration,patchLong,encounterLabel,...
    'VariableNames',{'worm #','encounter #','time enter (min)',...
    'duration (min)','log10(duration)','p(long duration)','label'});
sourceTableE = sortrows(sourceTableE,{'worm #','encounter #'});
writeSourceData([saveDir,'fig1_data2.xlsx'],{'Figure 1 - Source Data 2';'Figure 1E - Duration of patch encounters'},...
    {'Summary Statistics (Gaussian Mixture Model)';'Histogram Counts';'Source Data'},...
    {summaryTableE;histogramTableE;sourceTableE},{summaryNoteE;[];[]})

% timeEnter = encounter.timeEnter(indEncounter & indUncensored)./60;
timeEnter = encounter.timeEnter(indEncounter)./60;
pTimeEnter = ranksum(timeEnter(patchLong < 0.5),timeEnter(patchLong >= 0.5),'tail','left')
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
exportgraphics(gca,[saveDir,'fig1f.pdf'])

% Source Data
summaryTableF = table([numWorms;numWorms],[sum(patchLong < 0.5);sum(patchLong >= 0.5)],...
    splitapply(@mean,timeEnter,1+(patchLong >= 0.5)),splitapply(@std,timeEnter,1+(patchLong >= 0.5)),...
    splitapply(@median,timeEnter,1+(patchLong >= 0.5)),splitapply(@(x) prctile(x,2.5),timeEnter,1+(patchLong >= 0.5)),...
    splitapply(@(x) prctile(x,97.5),timeEnter,1+(patchLong >= 0.5)),{pTimeEnter;[]},...
    'VariableNames',{'# worms','# encounters','mean','std','median','2.5%','97.5%','p value'},'RowNames',{'short','long'},...
    'DimensionNames',{'encounter type','Variables'});
summaryNoteF = {'* mean, std, median, and 95% CI are given in minutes';'* p value calculated using one-tailed Mann-Whitney U-Test'};
writeSourceData([saveDir,'fig1_data3.xlsx'],{'Figure 1 - Source Data 3';'Figure 1F - Start time of patch encounters'},...
    {'Summary Statistics';'Source Data'},{summaryTableF;sourceTableE},{summaryNoteF;[]})

%% Figure 1G - Probability of Residing on Patch over Time (OD600 = 10)

oneHour = (0:3599)./60;
nReps = size(onPatchShuffled,3);
smoothSize = 180; % seconds

% Calculate on-patch (midpoint w/in 0.5 mm) 
onPatchUnshuffled = false(numWorms,length(oneHour));
for i = 1:numWorms
    indData = data.wormNum == wormNums(i);
    indInfo = info.plateNum == unique(data.plateNum(indData));
    onPatchUnshuffled(i,:) = interp1(data.timeOffset(indData),...
        double(data.distanceLawnEdge(indData) >= thresh.distMidpointEnter.*info.scale(indInfo)),...
        oneHour.*60) >= 0.5;
end

% Bootstrap sample unpermuted data to obtain confidence intervals
onPatchUnshuffledBoot = sort(bootstrp(nReps,@(X) mean(X,'omitnan'),onPatchUnshuffled));
onPatchUnshuffledCI = onPatchUnshuffledBoot(round(nReps.*[0.025,0.5,0.975]),:);

% Get pseudorandom patch position data
onPatchShuffledExp = onPatchShuffled(ismember(unique(data.wormNum),wormNums),:,:);
onPatchShuffledMean = sort(reshape(mean(onPatchShuffledExp,1,'omitnan'),[],nReps)');
onPatchShuffledCI = onPatchShuffledMean(round(nReps.*[0.025,0.5,0.975]),:);

% Build contingency table
contingency = reshape([sum(onPatchUnshuffled == 0,1);sum(onPatchShuffledExp == 0,[1,3]);...
    sum(onPatchUnshuffled == 1,1);sum(onPatchShuffledExp == 1,[1,3])],2,2,[]); % rows = {unshuffled;shuffled} , cols = {no,yes}
pPatchNoOverlap = nan(size(oneHour));
for i = 1:length(oneHour)
    % Fisher's Exact Test
    [~,pPatchNoOverlap(i)] = fishertest(contingency(:,:,i),'Tail','left');
    
    % Mann-Whitney U-Test
    % pPatchNoOverlap(i) = ranksum(onPatchUnshuffled(:,i),reshape(onPatchShuffled(:,i,:),[],1));
end

% Benjamini-Hochberg Correction
[hPatchNoOverlapBH,pPatchNoOverlapBH] = benjaminiHochberg(pPatchNoOverlap,0.001);
hPatchNoOverlapBH = double(hPatchNoOverlapBH);
hPatchNoOverlapBH(hPatchNoOverlapBH == 0) = NaN;

% Calculate average shuffled patch encounter
chancePatchEncounter = mean(onPatchShuffledCI(2,:))

figure; hold on;
patch([oneHour,flip(oneHour)],[movmean(onPatchShuffledCI(1,:),smoothSize,2),...
    movmean(flip(onPatchShuffledCI(3,:)),smoothSize,2)],...
    min(permuteColor + (1 - permuteColor).*(colorValue-0.1),1),'EdgeColor','none');
plot(oneHour,movmean(onPatchShuffledCI(2,:),smoothSize,2),'Color',permuteColor,'LineWidth',0.75);
patch([oneHour,flip(oneHour)],[movmean(onPatchUnshuffledCI(1,:),smoothSize,2),...
    movmean(flip(onPatchUnshuffledCI(3,:)),smoothSize,2)],'k','EdgeColor','none',...
    'FaceAlpha',1-colorValue);
plot(oneHour,movmean(onPatchUnshuffledCI(2,:),smoothSize,2),'k','LineWidth',0.75);
text(mean(oneHour(~isnan(hPatchNoOverlapBH))),0.97,'***',...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
plot(oneHour,hPatchNoOverlapBH.*1.05,'k','LineWidth',1)
xlim([0 60]); ylim([0 1.2]); xticks(0:15:60); yticks([0 1])
xlabel('time (min)','FontSize',8); ylabel('p(on patch)','FontSize',8)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('time (min)','FontSize',8); ylabel('p(on patch)','FontSize',8)
exportgraphics(gca,[saveDir,'fig1g.pdf'])

% Source data
summaryTableG = table([numWorms;numWorms],[nReps;nReps],{'bootstrapped samples';'permutations of patch location'},...
    'VariableNames',{'# worms','# replicates','replicate type'},'RowNames',{'original','permuted'},...
    'DimensionNames',{'patch locations','Variables'});
summaryNoteG = {'* patch locations were semi-randomly permuted as described in Figure 1 - supplement 4'};
contingencyTableG = array2table([oneHour',reshape(permute(contingency,[2,1,3]),4,[])',pPatchNoOverlap',pPatchNoOverlapBH'],...
    'VariableNames',{'time (min)','original - off patch','original - on patch','permuted - off patch','permuted - on patch','p-value','adjusted p-value'});
contingencyNoteG = {['* animal "on patch" if midpoint distance was within ',num2str(-thresh.distMidpointEnter),' mm from patch edge)'];...
    '* p-value computed using one-tailed Fisher''s Exact Test and adjusted using Benjamini-Hochberg method for multiple comparisons'};
sourceDataG = array2table([oneHour',onPatchUnshuffledCI([2,1,3],:)',onPatchShuffledCI([2,1,3],:)'],...
    'VariableNames',{'time (min)','median (original)','2.5% (original)','97.5% (original)',...
    'median (permuted)','2.5% (permuted)','97.5% (permuted)'});
writeSourceData([saveDir,'fig1_data4.xlsx'],{'Figure 1 - Source Data 4';'Figure 1G - On patch residence'},...
    {'Summary Statistics';'Contingency Tables (i.e. # replicates on or off patch given original and permuted patch locations)';...
    'Source Data (Median + 95% confidence interval for distribution of means of replicates)'},...
    {summaryTableG;contingencyTableG;sourceDataG},{summaryNoteG;contingencyNoteG;[]})

%% Diff of border amplitude (Growth vs. OD600 = 10)

borderDiff = mean(encounter.borderAmplitudeGrowth(indEncounter)./encounter.borderAmplitude(indEncounter))

%% Figure 1H - Example Worm Trace (Velocity)

wormNum = 360;
% plotTracks(info,data,wormNum,'velocitySmooth',saveDir,bodyPart,'yes');

%% Figure 1I,J,K -  Example Worm Velocity

wormNum = 360;

figure('Position',[360 278 840 420]);
hold on
for i = 1:length(indEncounterWorm)
    patch([encounter.timeEnter(indEncounterWorm(i)).*[1 1] ...
        encounter.timeExit(indEncounterWorm(i)).*[1 1]],...
        [0 400 400 0],colorValue.*[1 1 1],...
        'EdgeColor','none');
    if i == 3
        %patch(encounter.timeEnter(indEncounterWorm(i)) + [0 0 40 40],...
        patch([encounter.timeEnter(indEncounterWorm(i)).*[1 1], encounter.timeExit(indEncounterWorm(i)).*[1 1]],...
        [0 400 400 0],min(sampleColor + (1 - sampleColor).*(colorValue-0.1),1),...
        'EdgeColor','none');
    elseif i == 16
        %patch(encounter.timeEnter(indEncounterWorm(i)) + [0 0 40 40],...
        patch([encounter.timeEnter(indEncounterWorm(i)).*[1 1], encounter.timeExit(indEncounterWorm(i)).*[1 1]],...
        [0 400 400 0],min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1),...
        'EdgeColor','none');
    end
end

plot(data.timeOffset(indWorm),data.velocitySmooth(indWorm),'k','LineWidth',0.5)
xlim([0 3600]); xticks(0:600:3600); xticklabels(0:10:60); 
ylim([0 400]); yticks(0:100:400); 
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2.5 0.75])
xlabel('time (min)','FontSize',8); ylabel('velocity (μm/s)','FontSize',8)
exportgraphics(gca,[saveDir,'fig1i.pdf'])

xlim(encounter.timeEnter(indEncounterWorm(3)) + [-20 40])
xticks(encounter.timeEnter(indEncounterWorm(3)) + [-20:20:40])
xticklabels(-20:20:40); xlabel('time to patch (s)','FontSize',8)
set(gca,'Position',[ax.Position(1:2) 1 0.75])
exportgraphics(gca,[saveDir,'fig1j.pdf'])

xlim(encounter.timeEnter(indEncounterWorm(16)) + [-20 40])
xticks(encounter.timeEnter(indEncounterWorm(16)) + [-20:20:40])
exportgraphics(gca,[saveDir,'fig1k.pdf'])

%% Figure 1L - Slowing at Patch Edge (Time)

indTrajectory = find(ismember(trajectory.enterInfo,...
    encounter(indEncounter,trajectory.enterInfo.Properties.VariableNames)));

timeWindow = -60:0.25:240; % s
trajectoryMean = mean(trajectory.enterVelocity(indTrajectory,:),'omitnan');
[~,timeMax] = max(trajectoryMean); timeMax = timeWindow(timeMax);
[~,timeMin] = min(trajectoryMean); timeMin = timeWindow(timeMin);
indDecelerate = find(timeWindow >= timeMax & timeWindow <= timeMin);

deceleration = nan(length(indTrajectory),1);
dX = timeWindow(indDecelerate) - mean(timeWindow(indDecelerate));
for i = 1:length(indTrajectory)
    dY = trajectory.enterVelocity(indTrajectory(i),indDecelerate) - ...
        mean(trajectory.enterVelocity(indTrajectory(i),indDecelerate),'omitnan');
    deceleration(i) = sum(dX.*dY,'omitnan')./sum(dX.^2);
end

nReps = 1000;
trajectoryShiftBoot = nan(length(indTrajectory),length(timeWindow),nReps);
decelerationShiftBoot = nan(length(indTrajectory),nReps);
for i = 1:length(indTrajectory)
    indShift = find(data.wormNum == trajectory.enterInfo.wormNum(indTrajectory(i)));
    indEnter = randperm(length(indShift),nReps) + min(indShift) - 1;
    for j = 1:nReps
        trajectoryShiftBoot(i,:,j) = interp1(data.timeOffset(indShift),...
            data.velocitySmooth(indShift),data.timeOffset(indEnter(j)) + timeWindow);
        dY = trajectoryShiftBoot(i,indDecelerate,j) - ...
            mean(trajectoryShiftBoot(i,indDecelerate,j),'omitnan');
        decelerationShiftBoot(i,j) = sum(dX.*dY,'omitnan')./sum(dX.^2);
    end
end
trajectoryShiftCI95 = prctile(trajectoryShiftBoot,[2.5 50 97.5],[1,3]);
decelerationShiftCI95 = prctile(decelerationShiftBoot,[2.5 50 97.5],'all');

indTrajectory = ismember(trajectory.enterInfo,...
    encounter(indEncounter,trajectory.enterInfo.Properties.VariableNames));
trajectoryWorm = cell2mat(arrayfun(@(w) mean(trajectory.enterVelocity(indTrajectory & ...
    trajectory.enterInfo.wormNum == w,:),1,'omitnan'),wormNums,'UniformOutput',false));
trajectoryWormShift = cell2mat(arrayfun(@(w) mean(trajectoryShiftBoot(...
    trajectory.enterInfo.wormNum(indTrajectory) == w,:,randperm(nReps,1)),1,'omitnan'),wormNums,'UniformOutput',false));

%%
figure; hold on;
plot(timeWindow,trajectoryWormShift,'Color',min(permuteColor + (1 - permuteColor).*(colorValue-0.1),1))
plot(timeWindow,mean(trajectoryWormShift,'omitnan'),'Color',permuteColor,'LineWidth',1)
xline(0,'k')
plot(timeWindow,trajectoryWorm,'Color',colorValue.*[1 1 1])
plot(timeWindow,mean(trajectoryWorm,'omitnan'),'k','LineWidth',1)
xlim([-20 40]); ylim([0 400]); xticks(-20:20:40); yticks(0:100:400); 
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
xlabel('time to patch (s)','FontSize',8); ylabel('velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975 0.75])
exportgraphics(gca,[saveDir,'fig1l.pdf'],'ContentType','vector')

% Source data
summaryTableL = table([numWorms;numWorms],[length(indEncounter);length(indEncounter)],{'N/A';nReps},...
    'VariableNames',{'# worms','# encounters','# replicates'},'RowNames',{'original','permuted'},...
    'DimensionNames',{'encounter alignment','Variables'});
summaryNoteL = {'* original: time windows aligned to the start of each patch encounter';...
    '* permuted: time windows aligned to a random time point'};
sourceDataL = table(timeWindow',mean(trajectoryWorm,'omitnan')',mean(trajectoryWormShift,'omitnan')',...
    'VariableNames',{'time to patch entry (s)','original','permuted'});
sourceDataL(sourceDataL{:,1} < -20 | sourceDataL{:,1} > 40,:) = [];
sourceNoteL = {'* original: mean of mean trajectory for each animal (black line)';
    '* permuted: mean of means for only one replicate of permutation is shown (dark pink line)'};
writeSourceData([saveDir,'fig1_data5.xlsx'],{'Figure 1 - Source Data 5';'Figure 1L - Velocity of animals upon encounter with patch edge'},...
    {'Summary Statistics';'Source Data'},{summaryTableL;sourceDataL},{summaryNoteL;sourceNoteL})

%% Figure 1M - Deceleration at patch edge

% pAll = sum(mean(decelerationShiftBoot,1) < mean(deceleration))./size(decelerationShiftBoot,2)
% pShort = sum(mean(decelerationShiftBoot(patchLong < 0.5),1) < mean(deceleration(patchLong < 0.5)))./size(decelerationShiftBoot,2)
% pLong = sum(mean(decelerationShiftBoot(patchLong >= 0.5),1) < mean(deceleration(patchLong >= 0.5)))./size(decelerationShiftBoot,2)
pAll = 3*ranksum(deceleration,decelerationShiftBoot(:),'tail','left');
pShort = 3*ranksum(deceleration(patchLong < 0.5),decelerationShiftBoot(:),'tail','left');
pLong = 3*ranksum(deceleration(patchLong >= 0.5),decelerationShiftBoot(:),'tail','left');

figure
subplot(211); hold on
patch([0 4 4 0],prctile(mean(decelerationShiftBoot,1),[2.5 2.5 97.5 97.5]),...
    min(permuteColor + (1 - permuteColor).*(colorValue-0.1),1),'EdgeColor','none');
yline(prctile(mean(decelerationShiftBoot,1),50),'Color',permuteColor,'LineWidth',1)
v = Violin({deceleration},1,'ViolinColor',{colorValue.*[1 1 1]},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
v = Violin({deceleration(patchLong < 0.5)},2,'ViolinColor',{sampleColor},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
v = Violin({deceleration(patchLong >= 0.5)},3,'ViolinColor',{exploitColor},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
text(1,17,'***','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
text(2,17,'***','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
text(3,17,'***','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
set(gca,'Units','inches','XTick',1:3,'XTickLabel',{'all','short','long'},...
    'FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975*45/65 0.75])
xlim([0.4 3.6]); ylim([-50 35])
ylabel('Δ velocity (μm/s^2)','FontSize',8)
exportgraphics(gca,[saveDir,'fig1m.pdf'])

% Source data
summaryTableM = table([numWorms;numWorms;numWorms;numWorms],...
    [length(indEncounter);sum(patchLong < 0.5);sum(patchLong >= 0.5);size(decelerationShiftBoot,1)],...
    {'N/A';'N/A';'N/A';nReps},...
    [mean(deceleration);mean(deceleration(patchLong < 0.5));mean(deceleration(patchLong >= 0.5));mean(decelerationShiftBoot,'all')],...
    [std(deceleration);std(deceleration(patchLong < 0.5));std(deceleration(patchLong >= 0.5));std(decelerationShiftBoot,[],'all')],...
    [median(deceleration);median(deceleration(patchLong < 0.5));median(deceleration(patchLong >= 0.5));median(decelerationShiftBoot,'all')],...
    [prctile(deceleration,2.5);prctile(deceleration(patchLong < 0.5),2.5);prctile(deceleration(patchLong >= 0.5),2.5);prctile(decelerationShiftBoot,2.5,'all')],...
    [prctile(deceleration,97.5);prctile(deceleration(patchLong < 0.5),97.5);prctile(deceleration(patchLong >= 0.5),97.5);prctile(decelerationShiftBoot,97.5,'all')],...
    {pAll;pShort;pLong;[]},...
    'VariableNames',{'# worms','# encounters','# replicates','mean','std','median','2.5%','97.5%','p value'},'RowNames',{'all','short','long','permuted'},...
    'DimensionNames',{'encounter type','Variables'});
summaryNoteM = {'* deceleration measured in μm/s^2 as described in methods and Figure 1 - supplement 5';...
    '* p value calculated using one-tailed Mann-Whitney U-Test with Bonferonni Correction for multiple comparisons'};
sourceTableM = table(wormList,encounterList,deceleration,decelerationShiftBoot(:,1),patchLong,encounterLabel,...
    'VariableNames',{'worm #','encounter #','Δ velocity (μm/s^2)','Δ velocity (permuted)','p(long duration)','label'});
sourceTableM = sortrows(sourceTableM,{'worm #','encounter #'})
sourceNoteM = {'* example deceleration values shown for one replicate (of 1000) of permuting temporal alignment of patch encounter'};
writeSourceData([saveDir,'fig1_data6.xlsx'],{'Figure 1 - Source Data 6';'Figure 1M - Deceleration of animals upon encounter with patch edge'},...
    {'Summary Statistics';'Source Data'},{summaryTableM;sourceTableM},{summaryNoteM;sourceNoteM})

%% Figure 1N,O - Sustained slowing on patch (Distance)

indDistance = ismember(distance.info.wormNum,wormNums);
distanceWindow = -14:0.05:4.5; % mm

figure; hold on;
plot(distanceWindow,distance.velocitySmooth(indDistance,:),'Color',colorValue.*[1 1 1])
plot(distanceWindow,mean(distance.velocitySmooth(indDistance,:),'omitnan'),'k','LineWidth',1)
xline(0,'k--');
xlim([-2.5 0.75]); ylim([0 400]); yticks(0:100:400); xticks([-2:0,0.75])
plot([0 0.75],[320 320],'k','LineWidth',1); plot([thresh.distMidpointEnter,-2.5],[320 320],'k','LineWidth',1);
text(mean([0,0.75]),325,'on','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
text(mean([thresh.distMidpointEnter,-2.5]),320,'off','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
xlabel('distance to patch (mm)','FontSize',8); ylabel('avg. velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
exportgraphics(gca,[saveDir,'fig1n.pdf'])

% Source data
summaryTableN = table([numWorms],[length(indEncounter)],...
    'VariableNames',{'# worms','# encounters'});
sourceDataN = table(distanceWindow',mean(distance.velocitySmooth(indDistance,:),'omitnan')',...
    std(distance.velocitySmooth(indDistance,:),'omitnan')',median(distance.velocitySmooth(indDistance,:),'omitnan')',...
    prctile(distance.velocitySmooth(indDistance,:),2.5)', prctile(distance.velocitySmooth(indDistance,:),97.5)',...
    'VariableNames',{'distance to patch edge (mm)','mean','std','median','2.5%','97.5%'});
sourceDataN(sourceDataN{:,1} < -2.5 | sourceDataN{:,1} > 0.75,:) = [];
sourceNoteN = {'* distance between the animal’s midbody position and the patch edge (positive = inside patch; negative = outside patch)';
    '* mean, std, median, and 95% CI of mean velocities (μm/s) of each animal computed every 50 μm'};
writeSourceData([saveDir,'fig1_data7.xlsx'],{'Figure 1 - Source Data 7';'Figure 1N - Velocity of animals upon encounter with patch edge'},...
    {'Summary Statistics';'Source Data'},{summaryTableN;sourceDataN},{[];sourceNoteN})
%%
velocityOff = distance.velocitySmooth(indDistance,distanceWindow <= thresh.distMidpointEnter);
residenceOff = distance.probReside(indDistance,distanceWindow <= thresh.distMidpointEnter);
velocityOn = distance.velocitySmooth(indDistance,distanceWindow >= 0);
residenceOn = distance.probReside(indDistance,distanceWindow >= 0);
velocityOffWorm = sum(velocityOff.*residenceOff./sum(residenceOff,2),2,'omitnan'); % weighted sum is equivalent to mean of raw data
velocityOnWorm = sum(velocityOn.*residenceOn./sum(residenceOn,2),2,'omitnan');
[~,pSlowOn] = ttest(velocityOffWorm,velocityOnWorm,'tail','right')

figure; hold on
v = Violin({velocityOffWorm},1,'ViolinColor',{colorValue.*[1 1 1]},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
v = Violin({velocityOnWorm},2,'ViolinColor',{colorValue.*[1 1 1]},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
set(gca,'Units','inches','XTick',1:2,'XTickLabel',{'off','on'},...
    'FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975*35/65 0.75])
text(1.5,300,'***','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',8,'FontName','Arial')
plot([1 2],[320 320],'k','LineWidth',1)
xlim([0.4 2.6]); ylim([0 400]); yticks(0:100:400)
ylabel('avg. velocity (μm/s)','FontSize',8)
exportgraphics(gca,[saveDir,'fig1o.pdf'])

% Source Data
summaryTableO = table([numWorms;numWorms],[mean(velocityOffWorm);mean(velocityOnWorm)],...
    [std(velocityOffWorm);std(velocityOnWorm)],[median(velocityOffWorm);median(velocityOnWorm)],...
    [prctile(velocityOffWorm,2.5);prctile(velocityOnWorm,2.5)],[prctile(velocityOffWorm,97.5);prctile(velocityOnWorm,97.5)],...
    {[];pSlowOn},'VariableNames',{'# worms','mean','std','median','2.5%','97.5%','p value'},...
    'RowNames',{'off','on'},'DimensionNames',{'patch','Variables'})
summaryNoteO = {'* mean, std, median, and 95% CI are given in μm/s';'* p value calculated using one-tailed paired-sample t-test'};
sourceTableO = table((1:numWorms)',velocityOffWorm,velocityOnWorm,sum(residenceOff,2),sum(residenceOn,2),'VariableNames',...
    {'worm #','mean velocity off patch (μm/s)','mean velocity on patch (μm/s)','p(off patch)','p(on patch)'});
sourceNoteO = {'* mean velocity computed for each worm across all time points';...
    '* probabilites of residing on and off patch are given for each worm'};
writeSourceData([saveDir,'fig1_data8.xlsx'],{'Figure 1 - Source Data 8';'Figure 1O - Velocity on and off patch'},...
    {'Summary Statistics';'Source Data'},{summaryTableO;sourceTableO},{summaryNoteO;sourceNoteO})
%%

%% Figure 1N,O - Change in Speed at Patch Edge

figure; 
subplot(211); hold on;
xline(0,'k--')
Violin({deceleration},1,'ViolinColor',{colorValue.*[1 1 1]},...
    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
xlim([1/6 11/6]); ylim([-50 20])
set(gca,'Units','inches','XTick',[],'FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.3 0.75])
ylabel('Δ velocity (μm/s^2)','FontSize',8);
%exportgraphics(gca,[saveDir,'fig1n.pdf'])

%%
subplot(212); hold on
patch([-1.5 1.8 1.8 -1.5],decelerationShiftCI95([1 1 3 3]),...
    min(permuteColor + (1 - permuteColor).*(colorValue-0.1),1),'EdgeColor','none');
yline(decelerationShiftCI95(2),'Color',permuteColor,'LineWidth',1)
gscatter(log10(encounter.duration(indEncounter)./60),deceleration,...
   patchLong >= 0.5,[sampleColor;exploitColor],'.',5);
%scatter(log10(encounter.duration(indEncounter)./60),deceleration,5,colorValue.*[1 1 1],'filled');
xlim([-1.5 1.8]); ylim([-50 20]); legend('off')
xticks(log10([2 10 60 600 3600]./60)); xticklabels({'1/30','1/6','1','10','60'});
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlabel('duration (min)','FontSize',8); ylabel('Δ velocity (μm/s^2)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
exportgraphics(gca,[saveDir,'fig1o.pdf'])

pSlowdownLong = sum(deceleration(patchLong >= 0.5) < decelerationShiftCI95(1))./sum(patchLong >= 0.5)
pSlowdownShort = sum(deceleration(patchLong < 0.5) < decelerationShiftCI95(1))./sum(patchLong < 0.5)


%%
% ind = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
%     ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
%     strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn') & ...
%     encounter.timeEnter > 0);
% timePoints = [15 30 45 60];
% encounterSlowdown = nan(numWorms,length(timePoints));
% encounterSlowdownAll = arrayfun(@(w) mean(encounter.velocityBeforeEnter(ind & encounter.wormNum == w) - ...
%     encounter.velocityAfterEnter(ind & encounter.wormNum == w),'omitnan'),wormNums);
% 
% figure; hold on
% bar(1,mean(encounterSlowdownAll,'omitnan'),'FaceColor',colorValue.*[1 1 1])
% scatter(0.75 + rand(numWorms,1)./2,encounterSlowdownAll,5,'k','filled')
% xlim([1/6 11/6]); ylim([0 400]); yticks(0:100:400); % yticklabels(0:-100:-400)
% set(gca,'Units','inches','XTick',[],'FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
% ax = gca; set(gca,'Position',[ax.Position(1:2) 0.3 0.75])
% ylabel('Δ velocity (μm/s)','FontSize',8);
% exportgraphics(gca,[saveDir,'fig1l.pdf'])

for i = 1:numWorms
    for j = 1:length(timePoints)
        indSlow = ind & encounter.timeEnter./60 >= timePoints(j) - 15 & ...
            encounter.timeEnter./60 < timePoints(j) & encounter.wormNum == wormNums(i);
        encounterSlowdown(i,j) = mean(encounter.velocityBeforeEnter(indSlow) - ...
            encounter.velocityAfterEnter(indSlow),'omitnan');
    end
end
    
encounterSlowdownMean = mean(encounterSlowdown,'omitnan');
encounterSlowdownSEM = std(encounterSlowdown,'omitnan')./sqrt(numWorms);
figure; hold on;
bar(timePoints,encounterSlowdownMean,'FaceColor',colorValue.*[1 1 1])
scatter(repmat(timePoints,numWorms,1) + 6*rand(numWorms,length(timePoints)) - 3,...
    encounterSlowdown,5,'k','filled')
plot(timePoints.*[1;1],encounterSlowdownMean + encounterSlowdownSEM.*[1;-1],'k',...
    'LineWidth',0.75)
xlim([5 70]); xticks(timePoints);
ylim([0 400]); yticks(0:100:400); %yticklabels(0:-100:-400)
xlabel('time (min)','FontSize',8); ylabel('Δ velocity (μm/s)','FontSize',8)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975 0.75])
exportgraphics(gca,[saveDir,'fig1m.pdf'])
%% Figure 1Mv2 - Change in Speed at Patch Edge over time
ind = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn')) & ...
    encounter.timeEnter > 0;

figure; hold on;
scatter(encounter.timeEnter(ind)./60,encounter.velocityBeforeEnter(ind) - ...
            encounter.velocityAfterEnter(ind),5,colorValue.*[1 1 1],'filled');
scatter(encounter.timeEnter(indEncounter(3))./60,encounter.velocityBeforeEnter(indEncounter(3)) - ...
            encounter.velocityAfterEnter(indEncounter(3)),5,...
            min(exploreColor + (1 - exploreColor).*(colorValue-0.1),1),'filled');
scatter(encounter.timeEnter(indEncounter(16))./60,encounter.velocityBeforeEnter(indEncounter(16)) - ...
            encounter.velocityAfterEnter(indEncounter(16)),5,...
            min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1),'filled');
mdl = fitlm(encounter.timeEnter(ind)./60,encounter.velocityBeforeEnter(ind) - ...
            encounter.velocityAfterEnter(ind));
plot(0:60,predict(mdl,(0:60)'),'k','LineWidth',1)
xlim([0 60]); xticks(0:15:60); ylim([0 400]); yticks(0:100:400)
xlabel('time (min)','FontSize',8); ylabel('Δ velocity (μm/s)','FontSize',8)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
exportgraphics(gca,[saveDir,'fig1m2.pdf'])

%% Figure 1F - Encounter Cumulative Duration Over Time (OD600 = 10)

numSec = 3601;
patchDurationWorm = nan(numWorms,numSec);
for i = 1:numWorms
    ind = find(encounter.wormNum == wormNums(i) & strcmp(encounter.expName,expName) & ...
        ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
        strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn')));
    for j = 1:length(ind)
        timeWindow = round(encounter.timeEnter(ind(j))):round(encounter.timeExit(ind(j)));
        patchDurationWorm(i,timeWindow+1) = encounter.duration(ind(j));
        %patchDuration(i,timeWindow+1) = 1:length(timeWindow);
    end
end
% figure;
% subplot(211); imagesc(patchDurationWorm)
% subplot(212); plot(0:3600,mean(patchDurationWorm,'omitnan')./60)

oneHour = (0:3600)./60;
figure; 
subplot(211); hold on;

nReps = 1000;
beta = nan(nReps,1);
patchDurationShiftBoot = nan(nReps,length(oneHour));
for i = 1:nReps
    patchDurationShift = cell2mat(cellfun(@(dur) circshift(dur,randi(numSec),2),...
        mat2cell(patchDurationWorm,ones(numWorms,1)),'UniformOutput',false));
    % plot(oneHour,mean(patchDurationShift,'omitnan')./60)
    mdl = fitlm(oneHour(oneHour <= 30),mean(patchDurationShift(:,oneHour <= 30),'omitnan')./60,'linear');
    beta(i) = mdl.Coefficients.Estimate('x1');
    patchDurationShiftBoot(i,:) = mean(patchDurationShift,'omitnan')./60;
end
patchDurationShiftBoot = sort(patchDurationShiftBoot);
patchDurationShiftMedian = movmean(patchDurationShiftBoot(round(nReps.*0.5),:),300,2);
patchDurationShiftCI95 = movmean(patchDurationShiftBoot(round(nReps.*[0.025,0.975]),:),300,2);

patch([oneHour,flip(oneHour)],[patchDurationShiftCI95(1,:),flip(patchDurationShiftCI95(2,:))],...
    min(permuteColor + (1 - permuteColor).*(colorValue-0.1),1),'EdgeColor','none');
plot(oneHour,patchDurationShiftMedian,'Color',permuteColor,'LineWidth',1)

% patchDurationMean = movmean(mean(patchDurationWorm./60,'omitnan'),60);
% patchDurationSEM = movmean(std(patchDurationWorm./60,'omitnan')./sqrt(numWorms),60);
patchDurationBoot = sort(bootstrp(nReps,@(X) mean(X,'omitnan'),patchDurationWorm./60));
patchDurationMedian = movmean(patchDurationBoot(round(nReps.*0.5),:),300,2);
patchDurationCI95 = movmean(patchDurationBoot(round(nReps.*[0.025,0.975]),:),300,2);

noOverlapCI95 = double(patchDurationCI95(2,:) <= patchDurationShiftCI95(1,:) | ...
    patchDurationCI95(1,:) >= patchDurationShiftCI95(2,:));
noOverlapCI95(noOverlapCI95 == 0) = NaN;

% patch([oneHour,flip(oneHour)],[patchDurationMean - patchDurationSEM,flip(patchDurationMean + patchDurationSEM)],...
%     colorValue.*[1 1 1],'EdgeColor','none');
patch([oneHour,flip(oneHour)],[patchDurationCI95(1,:),flip(patchDurationCI95(2,:))],...
    colorValue.*[1 1 1],'EdgeColor','none');
plot(oneHour,patchDurationMedian,'k','LineWidth',1)
plot(oneHour,35*noOverlapCI95,'k','LineWidth',1)
maxNoOverlap = max(oneHour(~isnan(noOverlapCI95)))
text(maxNoOverlap/2,32,'*','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',8,'FontName','Arial')
mdl = fitlm(oneHour(oneHour <= 30),mean(patchDurationWorm(:,oneHour <= 30),'omitnan')./60,'linear');
% plot(oneHour(oneHour <= 30)',predict(mdl,oneHour(oneHour <= 30)'),'r','LineWidth',2)
xline(30,'k--')
xlim([0 60]); xticks(0:15:60); ylim([0 40])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975 0.75])
xlabel('time (min)','FontSize',8); ylabel('avg. duration (min)','FontSize',8)
exportgraphics(gca,[saveDir,'fig1f.pdf'])

subplot(212)
histogram(beta)
xline(mdl.Coefficients.Estimate('x1'),'r')
pStat = sum(beta > mdl.Coefficients.Estimate('x1'))./nReps

%%
% timePoints = [15 30 45 60];
% patchDurationMean = mean(patchDuration(:,timePoints.*60)./60,'omitnan');
% patchDurationSEM = std(patchDuration(:,timePoints.*60)./60,'omitnan')./sqrt(numWorms);
% subplot(4,4,[9,10,13,14]); hold on;
% bar(timePoints,patchDurationMean,'FaceColor',colorValue.*[1 1 1])
% scatter(repmat(timePoints,numWorms,1) + 6*rand(numWorms,length(timePoints)) - 3,...
%     patchDuration(:,timePoints.*60)./60,5,'k','filled')
% plot(timePoints.*[1;1],patchDurationMean + patchDurationSEM.*[1;-1],'k',...
%     'LineWidth',0.75)
% % b = boxchart(patchDuration(:,timePoints.*60)./60,'BoxFaceColor',colorValue.*[1 1 1],...
% %     'BoxFaceAlpha',1,'BoxEdgeColor','k');
% xlim([5 70]); xticks(timePoints); xlabel('time (min)','FontSize',8);
% ylabel('patch duration (min)','FontSize',8)
% set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
% ax = gca; set(gca,'Position',[ax.Position(1:2) 0.975 0.75])

%%
ind = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn'));

encounter.censorDuration = encounter.timeEnter == 0 | isnan(encounter.velocityOff); % right-censored (duration >= D)
encounter.censorEnter = -(encounter.timeEnter == 0); % left-censored (start <= T)
figure;
subplot(211); hold on;
[f,x] = ecdf(encounter.duration(ind)./60,'Censoring',encounter.censorDuration(ind),'Function','cdf','Bounds','on','Alpha',0.01);
ecdfhist(f,x)

subplot(212); hold on;
ecdf(encounter.timeEnter(ind)./60,'Censoring',encounter.censorEnter(ind),'Function','cdf','Bounds','on','Alpha',0.01)

figure;
scatter(encounter.timeEnter(ind)./60,encounter.timeExit(ind)./60)
%%
ind = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn'));
figure
scatter(encounter.timeEnter(ind)./60,log10(encounter.duration(ind)))

%%
% onPatchTime = repmat(oneHour,size(onPatchData,1),1);
% mdl = fitglm(onPatchTime(:),onPatchData(:),'Distribution','binomial');
% hold on; plot((0:60)',predict(mdl,(0:60)'))

chanceEncounter = nan(numWorms,1);
for i = 1:numWorms
    indInfo = info.plateNum == unique(data.plateNum(data.wormNum == wormNums(i)));
    patchArea = sum(bwdist(info.lawnMask{indInfo}) <= ...
        -info.scale(indInfo)*thresh.distMidpointEnter,'all');
    arenaArea = sum(info.arenaMask{indInfo},'all');
    chanceEncounter(i) = patchArea./arenaArea;
end
yline(mean(chanceEncounter),'k--','LineWidth',0.75)
text(0.95*max(oneHour),mean(chanceEncounter) + 0.05,...
    ['chance = ',num2str(100*mean(chanceEncounter),'%.1f'),'%'],...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'FontName','Arial','FontSize',8)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
exportgraphics(gca,[saveDir,'fig1g.pdf'])

p50 = oneHour(find(onPatchMean >= 0.5,1,'first'))
p95 = oneHour(find(onPatchMean >= 0.9,1,'first'))

%% Figure 1N - Encounters w/ Detectable Slowdown

% ind = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
%     ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
%     strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn'));
% 
% numSlowdowns = nan(numWorms,1);
% numEncounters = nan(numWorms,1);
% for i = 1:numWorms
%     numSlowdowns(i) = sum(encounter.wormNum == wormNums(i) & ind & ...
%         ~strcmp(encounter.label,'searchOn'));
%     numEncounters(i) = sum(encounter.wormNum == wormNums(i) & ind);
% end
% 
% figure; hold on
% bar(1,mean(numSlowdowns./numEncounters),'FaceColor',colorValue.*[1 1 1])
% scatter(0.75 + rand(numWorms,1)./2,numSlowdowns./numEncounters,5,'k','filled')
% plot([1;1],mean(numSlowdowns./numEncounters) + std(numSlowdowns./numEncounters)./sqrt(numWorms).*[1;-1],'k',...
%     'LineWidth',0.75)
% xlim([1/6 11/6])
% set(gca,'Units','inches','XTick',[],'FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
% ax = gca; set(gca,'Position',[ax.Position(1:2) 0.3 0.75])
% ylabel('p(slowing)','FontSize',8)
% exportgraphics(gca,[saveDir,'fig1n.pdf'])
% 
% slowdownCutOff = max(encounter.velocityBeforeEnter(strcmp(encounter.label,'searchOn')) - ...
%     encounter.velocityAfterEnter(strcmp(encounter.label,'searchOn')));
