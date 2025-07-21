%% Load data into workspace

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'labelEncounters.mat'),'exploitModel','senseModel',...
    'exploitDist','exploreDist');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
load(fullfile(path,expName,'trajectory.mat'),'trajectory');
load(fullfile(path,expName,'distance.mat'),'distance');
load(fullfile(path,expName,'permutePatchLocation.mat'),'onPatchShuffled');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
thresh = readtable([path,'foragingMini\encounterThresholds.csv']);
saveDir = [path,'figures\Fig2\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
searchOffColor = [240 228 66]./255; % yellow
permuteColor = [204 121 167]./255; % pink
exploreColor = [213 94 0]./255; % red

%% Get condition ids

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume','growthCondition',...
    'lawnOD600','peptone'}));
conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 0.5 & ...
    ~(GID.lawnOD600 == 0.1 & GID.growthCondition == 48));
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
conditionLabels = round([0;relativeBorder(2:end)],1,'significant');

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

%% Figure 2A - Concentrations tested w/ effective OD600

figure; hold on;
xValues = 1:numGroups; xValues(end-1:end) = xValues(end-1:end) + 0.5;
bar(xValues,relativeBorder,'FaceColor','flat','CData',colorValue.*[1 1 1])
for i = 1:numGroups
    t = text(xValues(i),relativeBorder(i),num2str(round(relativeBorder(i),1,'significant')),...
        'FontName','Arial','FontSize',6,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    if i == 1
        t.String = '0'; t.Position = [1 0.05 0];
    end
end
text(0.75,300,'Ø 1.7 mm','HorizontalAlignment','left','VerticalAlignment','middle',...
    'FontSize',8,'FontName','Arial')
xlim([0.25 numGroups+1.25])
xticks(xValues); xticklabels({'0','0.05','0.1','0.5','1','2',...
    '3','4','5','10','1','1'})
yscale('log'); ylim([0.05 1600]); yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 2 0.75])
ylabel('relative density','FontSize',8)
exportgraphics(gca,[saveDir,'fig2a.pdf'])

% Source data
relativeBorder_std = splitapply(@(X) std(X,'omitnan'),10*encounter.borderAmplitude./borderAmp10,G);
relativeBorder_std = relativeBorder_std(conditionG);
summaryTableA = table({'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',...
    [repmat(1,10,1);12;48],relativeBorder,relativeBorder_std,conditionLabels,colorValue,...
    'VariableNames',{'OD600','growth time (hr)','mean relative density','std relative density','label','saturation'})
summaryNoteA = {'* relative density was estimated as described in methods and Figure 2 - supplement 1-3';...
    '* bacteria-free patches were assigned a relative density 0.01 for use in log-based calculations';...
    '* condition labels use rounded relative density estimates for ease of reading';...
    '* saturation gives the grayscale color value from 0 (black) to 1 (white) for each condition'};
writeSourceData([saveDir,'fig2_data1.xlsx'],{'Figure 2 - Source Data 1';'Figure 2A - Relative density of bacterial patches in each condition'},...
    {'Summary Statistics'},{summaryTableA;},{summaryNoteA})

%% Figure 2B - Velocity for all OD as a function of dist. from patch edge

indDistance = ismember(distance.info.wormNum,wormNums);
distanceWindow = -14:0.05:4.5; % mm
distanceMean = splitapply(@(X) mean(X,'omitnan'),distance.velocitySmooth(indDistance,:),wormGroupInd);

figure; hold on;
plot(distanceWindow,distanceMean,'LineWidth',1)
colororder(colorValue.*[1 1 1])
xline(0,'k');
xlim([-2.5 0.75]); ylim([0 400]); yticks(0:100:400); xticks([-2:0,0.75])
plot([0 0.75],[320 320],'k','LineWidth',1); plot([thresh.distMidpointEnter,-2.5],[320 320],'k','LineWidth',1);
text(mean([0,0.75]),325,'on','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
text(mean([thresh.distMidpointEnter,-2.5]),320,'off','HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
xlabel('distance to patch (mm)','FontSize',8); ylabel('avg. velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
exportgraphics(gca,[saveDir,'fig2b.pdf'])

% Source data
summaryTableB = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',[ones(10,1);12;48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,...
    'VariableNames',{'label','OD600','growth time (hr)','# worms','saturation'})
summaryNoteB = {'* saturation gives the grayscale color value from 0 (black) to 1 (white) for each condition'};
sourceTableB = array2table([distanceWindow',distanceMean'],'VariableNames',...
    [{'distance to patch edge (mm)'},arrayfun(@num2str,conditionLabels,'UniformOutput',false)']);
sourceTableB(sourceTableB{:,1} < -2.5 | sourceTableB{:,1} > 0.75,:) = [];
sourceNoteB = {'* distance between the animal’s midbody position and the patch edge (positive = inside patch; negative = outside patch)';
    '* mean of mean velocities (μm/s) of each animal computed every 50 μm'};
writeSourceData([saveDir,'fig2_data2.xlsx'],{'Figure 2 - Source Data 2';'Figure 2B - Velocity of animals upon encounter with patch edge'},...
    {'Summary Statistics';'Source Data'},{summaryTableB;sourceTableB},{summaryNoteB;sourceNoteB})

%% Figure 2C - Velocity on patch for all worms

% Use residence to equally weigh data points to computer avg. velocity on patch
velocityOn = distance.velocitySmooth(indDistance,distanceWindow >= 0);
residenceOn = distance.probReside(indDistance,distanceWindow >= 0);
velocityOnWorm = sum(velocityOn.*residenceOn./sum(residenceOn,2),2,'omitnan');
velocityOnMean = splitapply(@(X) mean(X,'omitnan'),velocityOnWorm,wormGroupInd);

% Mann-Whitney U-Test w/ Bonferoni correction to compare velocity on patch
% to velocity on OD = 0 patches
velocityOnZero = velocityOnWorm(wormGroupInd == find(GID.lawnOD600(conditionG) == 1e-10));
pSlowOn = numGroups*splitapply(@(X) ranksum(X,velocityOnZero,'tail','left'),velocityOnWorm,wormGroupInd)

figure; hold on
gscatter(wormBorder,velocityOnWorm,wormGroup,colorValue.*[1 1 1],'.',5)
[f,g] = fit(log10(wormBorder)-log10(borderAmp0)+1,velocityOnWorm,'logistic4');
x = logspace(log10(borderAmp0)-1,log10(max(relativeBorder))+1,100);
plot(x,f(log10(x) - log10(borderAmp0) + 1),'k','LineWidth',1)
scatter(relativeBorder,velocityOnMean,20,colorValue.*[1 1 1],'filled','MarkerEdgeColor','k')
xscale('log'); xlim([0.005 1000]); ylim([0 400]); yticks(0:100:400)
xticks([0.01 0.1 1 10 100]); xticklabels([0 0.1 1 10 100]); legend('off')
text(10^mean(log10(prctile(relativeBorder(pSlowOn < 1e-3),[0 100]))),300,'***',...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',8,'FontName','Arial')
plot(prctile(relativeBorder(pSlowOn < 1e-3),[0 100]),[320 320],'k','LineWidth',1)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlabel('relative density','FontSize',8); ylabel('velocity on (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 0.975 0.75])
exportgraphics(gca,[saveDir,'fig2c.pdf'])

% Source data
summaryTableC = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',...
    [repmat(1,10,1);12;48],relativeBorder,...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,...
    velocityOnMean,[{[]};num2cell(pSlowOn(2:end))],...
    'VariableNames',{'label','OD600','growth time (hr)','relative density','# worms','saturation','mean velocity on patch (μm/s)','p value'});
summaryNoteC = {'* mean of mean velocity (μm/s) on patch for each animal';...
    '* p value calculated by comparing to 0 condition (i.e. bacteria-free patch) using one-tailed Mann Whitney U-Test with Bonferroni Correction for multiple comparisons'};
sigmoidTableC = table(f.a,f.b,f.c,f.d,g.rsquare,'VariableNames',...
    {'a','b','c','d','r-squared'});
sourceTableC = table((1:numWorms)',conditionLabels(wormGroupInd),velocityOnWorm,sum(residenceOn,2),'VariableNames',...
    {'worm #','label','mean velocity on patch (μm/s)','p(on patch)'});
sourceNoteC = {'* mean velocity computed for each worm across all time points';...
    '* probability of residing on patch is given for each worm'};
writeSourceData([saveDir,'fig2_data3.xlsx'],{'Figure 2 - Source Data 3';'Figure 2C - Velocity on patch'},...
    {'Summary Statistics';'Sigmoid (y = d + (a-d)/(1 + (x/c)^b)';'Source Data'},{summaryTableC;sigmoidTableC;sourceTableC},{summaryNoteC;[];sourceNoteC})

%% Figure 2D - Example traces OD600 = 0, 1, 10, 1(48H)

wormNum = [465, 491, 343, 360, 624]; % [0, 1, 4, 10, 1(48H)]
% for i = 1:length(wormNum)
%     plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
% end

%% Figure 2E - Time on Patch as a function of border amplitude

timeOnPatch = arrayfun(@(w) sum(encounter.duration(ismember(G,conditionG) & ...
    encounter.wormNum == w & indEncounter))/60,wormNums);
timeOnMean = splitapply(@(X) mean(X,'omitnan'),timeOnPatch,wormGroupInd);

figure; hold on
gscatter(wormBorder,timeOnPatch,wormGroup,colorValue.*[1 1 1],'.',5)
[f,g] = fit(log10(wormBorder)-log10(borderAmp0)+1,timeOnPatch,'logistic4','Upper',[Inf Inf Inf 60]);
x = logspace(log10(borderAmp0)-1,log10(max(relativeBorder))+1,100);
plot(x,f(log10(x) - log10(borderAmp0) + 1),'k','LineWidth',1)
scatter(relativeBorder,timeOnMean,20,colorValue.*[1 1 1],'filled','MarkerEdgeColor','k')
xscale('log'); xlim([0.005 1000]); ylim([0 60]); yticks(0:20:60)
xticks([0.01 0.1 1 10 100]); xticklabels([0 0.1 1 10 100]); legend('off')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlabel('relative density','FontSize',8); ylabel('time on patch (min)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 0.975 0.75])
exportgraphics(gca,[saveDir,'fig2e.pdf'])

% Kendall's correlation (tau is + for monotonically increasing y)
[tauTimeIncrease,pTimeIncrease] = corr(wormBorder,timeOnPatch,'type','Kendall')

% Source data
summaryTableE = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',...
    [repmat(1,10,1);12;48],relativeBorder,...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,...
    timeOnMean,[{tauTimeIncrease};repmat({[]},numGroups-1,1)],[{pTimeIncrease};repmat({[]},numGroups-1,1)],...
    'VariableNames',{'label','OD600','growth time (hr)','relative density','# worms','saturation','mean time on patch (min)','τ','p value'});
summaryNoteE = {'* mean of mean time on patch (min) for each animal';...
    '* p value calculated using Kendall''s τ rank correlation coefficient (test of monotonicity)'};
sigmoidTableE = table(f.a,f.b,f.c,f.d,g.rsquare,'VariableNames',...
    {'a','b','c','d','r-squared'});
sourceTableE = table((1:numWorms)',conditionLabels(wormGroupInd),timeOnPatch,...
    'VariableNames',{'worm #','label','time on patch (min)'});
sourceNoteE = {'* mean time on patch computed for each worm across all time points';...
    '* probability of residing on patch is given for each worm'};
writeSourceData([saveDir,'fig2_data4.xlsx'],{'Figure 2 - Source Data 4';'Figure 2E - Time on patch'},...
    {'Summary Statistics';'Sigmoid (y = d + (a-d)/(1 + (x/c)^b)';'Source Data'},{summaryTableE;sigmoidTableE;sourceTableE},{summaryNoteE;[];sourceNoteE})

%% Figure 2F - p(on patch) as function of time

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
onPatchUnshuffledBoot = nan(nReps,length(oneHour),numGroups);
for i = 1:numGroups
    onPatchUnshuffledBoot(:,:,i) = sort(bootstrp(nReps,@(X) mean(X,'omitnan'),onPatchUnshuffled(wormGroupInd == i,:)));
end
onPatchUnshuffledCI = onPatchUnshuffledBoot(round(nReps.*[0.025,0.5,0.975]),:,:);

% Get pseudorandom patch position data
onPatchShuffledExp = onPatchShuffled(ismember(unique(data.wormNum),wormNums),:,:);
onPatchShuffledMean = nan(nReps,length(oneHour),numGroups);
for i = 1:numGroups
     onPatchShuffledMean(:,:,i) = sort(reshape(mean(onPatchShuffledExp(wormGroupInd == i,:,:),1,'omitnan'),[],nReps)',1);
end
onPatchShuffledCI = onPatchShuffledMean(round(nReps.*[0.025,0.5,0.975]),:,:);

% Build contingency table
pPatchNoOverlap = nan(length(oneHour),numGroups);
pPatchNoOverlapBH = nan(length(oneHour),numGroups);
hPatchNoOverlapBH = nan(length(oneHour),numGroups);
contingency = nan(2,2,length(oneHour),numGroups);
for j = 1:numGroups
    contingency(:,:,:,j) = reshape([sum(onPatchUnshuffled(wormGroupInd == j,:) == 0,1);sum(onPatchShuffledExp(wormGroupInd == j,:,:) == 0,[1,3]);...
        sum(onPatchUnshuffled(wormGroupInd == j,:) == 1,1);sum(onPatchShuffledExp(wormGroupInd == j,:,:) == 1,[1,3])],2,2,[]); % rows = {unshuffled;shuffled} , cols = {no,yes}
    for i = 1:length(oneHour)
        % Fisher's Exact Test
        [~,pPatchNoOverlap(i,j)] = fishertest(contingency(:,:,i,j),'Tail','left');

        % Mann-Whitney U-Test
        % pPatchNoOverlap(i) = ranksum(onPatchUnshuffled(:,i),reshape(onPatchShuffled(:,i,:),[],1));
    end

    % Benjamini-Hochberg Correction
    [hPatchNoOverlapBH(:,j),pPatchNoOverlapBH(:,j)] = benjaminiHochberg(pPatchNoOverlap(:,j),0.001);
end
hPatchNoOverlapBH = double(hPatchNoOverlapBH);
hPatchNoOverlapBH(hPatchNoOverlapBH == 0) = NaN;

figure; hold on
%plot(oneHour,movmean(onPatchShuffledMeanMean,1,1),'LineWidth',0.75);
plot(oneHour,permute(movmean(onPatchUnshuffledCI(2,:,:),smoothSize,2),[2,3,1]),'LineWidth',1);
plot(oneHour'.*hPatchNoOverlapBH,permute(movmean(onPatchUnshuffledCI(2,:,:),smoothSize,2),[2,3,1]),'LineWidth',1);
colororder([min(permuteColor + (1 - permuteColor).*(colorValue-0.1),1);colorValue.*[1 1 1]])
xlim([0 60]); ylim([0 1]); xticks(0:15:60)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('time (min)','FontSize',8); ylabel('p(on patch)','FontSize',8)
exportgraphics(gca,[saveDir,'fig2f.pdf'])

% Source data
[~,earliestH] = max(hPatchNoOverlapBH,[],1);
earliestH = oneHour(earliestH);
earliestH(mean(isnan(hPatchNoOverlapBH),1) == 1) = NaN;
summaryTableF = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',...
    [repmat(1,10,1);12;48],relativeBorder,...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,earliestH',...
    'VariableNames',{'label','OD600','growth time (hr)','relative density','# worms','saturation','earliest time (min) p<0.001'});
summaryNoteF = {'* patch locations were semi-randomly permuted (1000 replicates) as described in Figure 1 - supplement 4';...
    '* the earliest time (min) when animals p(on patch) exceeded that of the permuted estimates is given'};
contingencyTableF = array2table([reshape(repmat(conditionLabels,1,length(oneHour))',[],1),...
    repmat(oneHour',numGroups,1),reshape(permute(contingency,[2,1,3,4]),4,[])',reshape(pPatchNoOverlap,1,[])',reshape(pPatchNoOverlapBH,1,[])'],...
    'VariableNames',{'label','time (min)','original - off patch','original - on patch','permuted - off patch','permuted - on patch','p-value','adjusted p-value'});
contingencyNoteF = {['* animal "on patch" if midpoint distance was within ',num2str(-thresh.distMidpointEnter),' mm from patch edge)'];...
    '* p-value computed using one-tailed Fisher''s Exact Test and adjusted using Benjamini-Hochberg method for multiple comparisons'};
sourceDataF = array2table([reshape(repmat(conditionLabels,1,length(oneHour))',[],1),...
    repmat(oneHour',numGroups,1),reshape(onPatchUnshuffledCI([2,1,3],:,:),3,[])',...
    reshape(onPatchShuffledCI([2,1,3],:,:),3,[])'],...
    'VariableNames',{'label','time (min)','median (original)','2.5% (original)','97.5% (original)',...
    'median (permuted)','2.5% (permuted)','97.5% (permuted)'});
writeSourceData([saveDir,'fig2_data5.xlsx'],{'Figure 2 - Source Data 5';'Figure 2F - On patch residence'},...
    {'Summary Statistics';'Contingency Tables (i.e. # replicates on or off patch given original and permuted patch locations)';...
    'Source Data (Median + 95% confidence interval for distribution of means of replicates)'},...
    {summaryTableF;contingencyTableF;sourceDataF},{summaryNoteF;contingencyNoteF;[]})

%% Figure 2G - Patch Duration

figure;  hold on
durationBins = -2:0.01:2;
kDense = nan(length(durationBins),numGroups);
for i = 1:numGroups
    kDense(:,i) = kde(log10(encounter.duration(indEncounter & G == conditionG(i))./60),...
        'EvaluationPoints',durationBins,'Bandwidth',4/20);
    plot(durationBins,kDense(:,i),'LineWidth',1,'Color',colorValue(i).*[1 1 1]);
end
xlim([-1.7 2]); ylim([0 1.7]); xticks(log10([2 10 60 600 3600]./60)); 
xticklabels({'1/30','1/6','1','10','60'});
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','YTick',[]); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('encounter duration (min)','FontSize',8); ylabel('p(encounter)','FontSize',8)
exportgraphics(gca,[saveDir,'fig2g.pdf'])

% Source data
summaryTableG = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',[ones(10,1);12;48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',arrayfun(@(g) sum(indEncounter & G == conditionG(g)),1:numGroups)',colorValue,...
    'VariableNames',{'label','OD600','growth time (hr)','# worms','# encounters','saturation'});
summaryNoteG = {'* saturation gives the grayscale color value from 0 (black) to 1 (white) for each condition'};
kdeTableG = array2table([10.^durationBins',durationBins',kDense],'VariableNames',...
    [{'duration (min)','log10(duration)'},arrayfun(@num2str,conditionLabels,'UniformOutput',false)']);
kdeNoteG = {'* kernel density estimates (KDE) of distributions of duration of patch encounter';...
    '* KDEs are reported in lieu of histograms for visibility purposes only'};
sourceTableG = arrayfun(@(i) histcounts(log10(encounter.duration(indEncounter & G == conditionG(i))./60),...
    -2.05:0.1:2.05)',1:numGroups,'UniformOutput',false);
sourceTableG = array2table([10.^(-2:0.1:2)',(-2:0.1:2)',horzcat(sourceTableG{:})],'VariableNames',...
    [{'duration (min)','log10(duration)'},arrayfun(@num2str,conditionLabels,'UniformOutput',false)']);
sourceNoteG = {'* # of encounters observed for each duration bin'};
writeSourceData([saveDir,'fig2_data6.xlsx'],{'Figure 2 - Source Data 6';'Figure 2G - Duration of patch encounters'},...
    {'Summary Statistics';'Kernel Density Estimate';'Source Data'},{summaryTableG;kdeTableG;sourceTableG},{summaryNoteG;kdeNoteG;sourceNoteG})

%% Figure 2H - Explore vs. Exploit Clusters

% Get contours of standard deviation(s) of gaussians
alphaVal = [0.4 0.3 0.2];
stdValues = [1 2 3];
numPoints = 1000;
exploitEllipses = gaussianContours(exploitDist.mu,exploitDist.Sigma,stdValues,numPoints);
exploreEllipses = gaussianContours(exploreDist.mu,exploreDist.Sigma,stdValues,numPoints);

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
for i = 1:length(alphaVal)
    patch(exploitEllipses(:,2*i-1),exploitEllipses(:,2*i),...
        min(exploitColor + (1 - exploitColor).*(colorValue(10)-0.1),1),...
        'FaceAlpha',alphaVal(i),'EdgeAlpha',0)
    patch(exploreEllipses(:,2*i-1),exploreEllipses(:,2*i),...
        min(sampleColor + (1 - sampleColor).*(colorValue(10)-0.1),1),...
        'FaceAlpha',alphaVal(i),'EdgeAlpha',0)
end
scatter(log10(encounter.duration(indEncounter)./60),log10(encounter.velocityOn(indEncounter)),1,'CData',...
    encounter.exploitPosterior(indEncounter).*exploitColor + ...
    (1-encounter.exploitPosterior(indEncounter)).*encounter.sensePosterior(indEncounter).*sampleColor + ...
    (1-encounter.exploitPosterior(indEncounter)).*(1-encounter.sensePosterior(indEncounter)).*searchOnColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none')

xlim([-2 2]); ylim([1.3 2.8]); legend('off')
xticks(log10([2 10 60 600 3600]./60)); xticklabels({'1/30','1/6','1','10','60'});
yticks(log10([25 50 100 200 400])); yticklabels([25 50 100 200 400])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 1.5 1.5])
xlabel('encounter duration (min)','FontSize',8); ylabel('avg. velocity on patch (μm/s)','FontSize',8)

labelNames = {'searchOn','sample','exploit'};
velocityBins = 1.3:0.01:2.8;
kDense = nan(length(labelNames),length(velocityBins));
kDense(1,:) = kde(log10(encounter.velocityOn(indEncounter)),...
        'EvaluationPoints',velocityBins,'Bandwidth',1.5/20,'Weight',...
        (1 - encounter.exploitPosterior(indEncounter)).*(1 - encounter.sensePosterior(indEncounter)));
kDense(2,:) = kde(log10(encounter.velocityOn(indEncounter)),...
        'EvaluationPoints',velocityBins,'Bandwidth',1.5/20,'Weight',...
        (1 - encounter.exploitPosterior(indEncounter)).*encounter.sensePosterior(indEncounter));
kDense(3,:) = kde(log10(encounter.velocityOn(indEncounter)),...
        'EvaluationPoints',velocityBins,'Bandwidth',1.5/20,'Weight',...
        encounter.exploitPosterior(indEncounter));
subplot(6,6,3 + [3 9 15]'); hold on
for i = 1:length(labelNames)
    plot(velocityBins,kDense(i,:),'LineWidth',1,'Color',eval([labelNames{i},'Color']));
end
set(gca,'View',[90 -90],'XLim',[1.3 2.8],'XTick',[],'YTick',[],...
    'YAxisLocation','right','Units','inches','FontName','Arial',...
    'FontSize',8,'Box','on','TickDir','out','Ylim',[0 1.1*max(kDense,[],'all')]);
ax = gca; set(gca,'Position',[4 4 0.5 1.5])

durationBins = -2:0.01:2;
kDense = nan(length(labelNames),length(durationBins));
kDense(1,:) = kde(log10(encounter.duration(indEncounter)./60),...
        'EvaluationPoints',durationBins,'Bandwidth',4/20,'Weight',...
        (1 - encounter.exploitPosterior(indEncounter)).*(1 - encounter.sensePosterior(indEncounter)));
kDense(2,:) = kde(log10(encounter.duration(indEncounter)./60),...
        'EvaluationPoints',durationBins,'Bandwidth',4/20,'Weight',...
        (1 - encounter.exploitPosterior(indEncounter)).*encounter.sensePosterior(indEncounter));
kDense(3,:) = kde(log10(encounter.duration(indEncounter)./60),...
        'EvaluationPoints',durationBins,'Bandwidth',4/20,'Weight',...
        encounter.exploitPosterior(indEncounter));
subplot(6,6,[13 14 15] + 18); hold on
for i = 1:length(labelNames)
    plot(durationBins,kDense(i,:),'LineWidth',1,'Color',eval([labelNames{i},'Color']));
end
set(gca,'XLim',[-2 2],'XTick',[],'YTick',[],'Units','inches','FontName','Arial',...
    'FontSize',8,'Box','on','TickDir','out','Ylim',[0 1.1*max(kDense,[],'all')]);
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 1.5 0.5])
exportgraphics(gcf,[saveDir,'fig2h.pdf'],'ContentType','vector')

% Source data
[~,wormList] = ismember(encounter.wormNum(indEncounter),wormNums);
[~,conditionList] = ismember(G(indEncounter),conditionG);
encounterList = arrayfun(@(w) 1:sum(encounter.wormNum(indEncounter)==w),unique(encounter.wormNum(indEncounter)),'UniformOutput',false);
encounterList = [encounterList{:}]';
summaryTableH = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',[ones(10,1);12;48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',arrayfun(@(g) sum(indEncounter & G == conditionG(g)),1:numGroups)',...
    'VariableNames',{'label','OD600','growth time (hr)','# worms','# encounters'});
gmmTableH = array2table([[exploreDist.mu;exploitDist.mu],...
    [exploreDist.Sigma([1,4,2,]);exploitDist.Sigma([1,4,2,])],...
    reshape(exploitModel.ComponentProportion,2,1)],...
    'VariableNames',{'mu - duration','mu - velocity','sigma - var(duration)',...
    'sigma - var(velocity)','sigma - cov','proportion'},'RowNames',{'explore','exploit'},...
    'DimensionNames',{'encounter type','Variables'});
gmmNoteH = {'* mu, sigma, and proportion define the mean, covariance matrix, and component proportion of each gaussian';...
    '* mu and sigma are in units of log10(encounter duration in minutes) and log10(velocity in μm/s)'};
sourceTableH = table(conditionLabels(conditionList),wormList,encounterList,...
    encounter.duration(indEncounter)./60,log10(encounter.duration(indEncounter)./60),...
    encounter.velocityOn(indEncounter),log10(encounter.velocityOn(indEncounter)),...
    encounter.exploitPosterior(indEncounter),1-encounter.exploitPosterior(indEncounter),...
    encounter.sensePosterior(indEncounter),1-encounter.sensePosterior(indEncounter),...
    'VariableNames',{'label','worm #','encounter #','encounter duration (min)',...
    'log10(duration)','mean velocity on patch (μm/s)','log10(velocity)',...
    'p(exploit)','p(explore)','p(sense)','p(non-sense)'});
sourceNoteH = {'* posterior probability of classification as explore or exploit as estimated by the Gaussian Mixture Model described in methods and Figure 2 - supplement 6';...
    '* posterior probability of classification as sense or non-sense as estimated by semi-supervised Quadratic Discriminant Analysis (QDA) described in methods, Figure 2I, Figure 2 - supplement 7, and Video 5';...
    '* RGB value calculation: p(exploit)*[0 114 178] + p(explore)*p(sense)*[0 158 115] + p(explore)*p(non-sense)*[230 159 0]'};
writeSourceData([saveDir,'fig2_data7.xlsx'],{'Figure 2 - Source Data 7';'Figure 2H - Encounter classification as explore or exploit using a Gaussian Mixture Model (GMM)'},...
    {'Summary Statistics';'Gaussian Mixture Model';'Source Data'},{summaryTableH;gmmTableH;sourceTableH},{[];gmmNoteH;sourceNoteH})

%% Figure 2I - Sample vs. Search Clusters

% Get QDA paraboloid boundary
K = senseModel.Coeffs(1,2).Const;
L = senseModel.Coeffs(1,2).Linear; 
Q = senseModel.Coeffs(1,2).Quadratic;
quadratic3D = @(x1,x2,x3) K + L(1)*x1 + L(2)*x2 + L(3)*x3 + ...
    Q(1,1)*x1.^2 + Q(2,2)*x2.^2 + Q(3,3)*x3.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + (Q(1,3)+Q(3,1))*x1.*x3 + (Q(2,3)+Q(3,2))*x2.*x3;
x2_median = prctile(encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),50);
quadratic2D_13_median = @(x1,x3) K + L(1)*x1 + L(2)*x2_median + L(3)*x3 + ...
    Q(1,1)*x1.^2 + Q(2,2)*x2_median.^2 + Q(3,3)*x3.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2_median + (Q(1,3)+Q(3,1))*x1.*x3 + (Q(2,3)+Q(3,2))*x2_median.*x3;

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
fp = fimplicit(quadratic2D_13_median,[0 450 -60 40],'LineStyle','none');
patch([450 fp.XData 450],[-60 fp.YData 40],min(searchOnColor + (1 - searchOnColor).*(colorValue(10)-0.1),1),...
         'FaceAlpha',min(alphaVal),'EdgeAlpha',0)
patch([0 fp.XData 0],[-60 fp.YData 40],min(sampleColor + (1 - sampleColor).*(colorValue(10)-0.1),1),...
         'FaceAlpha',min(alphaVal),'EdgeAlpha',0)
scatter(encounter.velocityOnMin(indEncounter),encounter.decelerate(indEncounter),...
    1,'CData',encounter.exploitPosterior(indEncounter).*exploitColor + ...
    (1-encounter.exploitPosterior(indEncounter)).*encounter.sensePosterior(indEncounter).*sampleColor + ...
    (1-encounter.exploitPosterior(indEncounter)).*(1-encounter.sensePosterior(indEncounter)).*searchOnColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none')
legend('off'); xlim([0 450]); ylim([-60 40]); 
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 1.5 1])
xlabel('min. velocity on patch (μm/s)','FontSize',8); ylabel('Δ velocity (μm/s^2)','FontSize',8)

labelNames = {'searchOn','sample','exploit'};
velocityBins = -60:40;
kDense = nan(length(labelNames),length(velocityBins));
kDense(1,:) = kde(encounter.decelerate(indEncounter),...
        'EvaluationPoints',velocityBins,'Bandwidth',100/20,'Weight',...
        (1 - encounter.exploitPosterior(indEncounter)).*(1 - encounter.sensePosterior(indEncounter)));
kDense(2,:) = kde(encounter.decelerate(indEncounter),...
        'EvaluationPoints',velocityBins,'Bandwidth',100/20,'Weight',...
        (1 - encounter.exploitPosterior(indEncounter)).*encounter.sensePosterior(indEncounter));
kDense(3,:) = kde(encounter.decelerate(indEncounter),...
        'EvaluationPoints',velocityBins,'Bandwidth',100/20,'Weight',...
        encounter.exploitPosterior(indEncounter));

subplot(6,6,3 + [3 9 15]'); hold on
for i = 1:length(labelNames)
    plot(velocityBins,kDense(i,:),'LineWidth',1,'Color',eval([labelNames{i},'Color']));
end
set(gca,'View',[90 -90],'XLim',[-60 40],'XTick',[],'YTick',[],...
    'YAxisLocation','right','Units','inches','FontName','Arial',...
    'FontSize',8,'Box','on','TickDir','out','Ylim',[0 1.1*max(kDense,[],'all')]);
ax = gca; set(gca,'Position',[4 4 0.5 1])
exportgraphics(gcf,[saveDir,'fig2i.pdf'],'ContentType','vector')

% Source data
summaryTableI = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',[ones(10,1);12;48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',arrayfun(@(g) sum(indEncounter & G == conditionG(g)),1:numGroups)',...
    'VariableNames',{'label','OD600','growth time (hr)','# worms','# encounters'});
qdaTableI = array2table([K,L',Q(1,1),Q(2,2),Q(3,3),Q(1,2)+Q(2,1),Q(1,3)+Q(3,1),Q(2,3)+Q(3,2)],...
    'VariableNames',{'a','b','c','d','e','f','g','h','i','j'});
qdaNoteI = {'* x1 = min. velocity on patch (μm/s); x2 = max. Δ velocity (μm/s^2); x3 = Δ velocity (μm/s^2)'};
sourceTableI = table(conditionLabels(conditionList),wormList,encounterList,...
    encounter.velocityOnMin(indEncounter),encounter.decelerate(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    encounter.exploitPosterior(indEncounter),1-encounter.exploitPosterior(indEncounter),...
    encounter.sensePosterior(indEncounter),1-encounter.sensePosterior(indEncounter),...
    rand(sum(indEncounter),1) <= encounter.exploitPosterior(indEncounter),conditionLabels(conditionList) == 0,...
    'VariableNames',{'label','worm #','encounter #','min. velocity on patch (μm/s)',...
    'Δ velocity (μm/s^2)','max. Δ velocity (μm/s^2)',...
    'p(exploit)','p(explore)','p(sense)','p(non-sense)',...
    'labeled sense (QDA)','labeled non-sense (QDA)'})
sourceNoteI = {'* posterior probability of classification as explore or exploit as estimated by the Gaussian Mixture Model described in methods, Figure 2H, and Figure 2 - supplement 6';...
    '* posterior probability of classification as sense or non-sense as estimated by semi-supervised quadratic discriminant analysis (QDA) described in methods, Figure 2 - supplement 7, and Video 5';...
    '* RGB value calculation: p(exploit)*[0 114 178] + p(explore)*p(sense)*[0 158 115] + p(explore)*p(non-sense)*[230 159 0]';...
    '* semi-supervied QDA uses labeled data to generate a better discriminator';...
    '    data labeled non-sense and sense are indicated';...
    '    data labeled sense were included probabilistically (e.g. if p(exploit) = 0.9, the encounter would be labeled sense for 90% of replicates); labels here represent one replicate (of 1000)'};
writeSourceData([saveDir,'fig2_data8.xlsx'],{'Figure 2 - Source Data 8';'Figure 2I - Encounter classification as sense or non-sense using semi-supervised Quadratice Discriminant Analysis (QDA)'},...
    {'Summary Statistics';'Quadratic Discriminant (y = a + b*x1 + c*x2 + d*x3 + e*x1^2 + f*x2^2 + g*x3^2 + h*x1*x2 + i*x1*x3 + j*x2*x3)';'Source Data'},...
    {summaryTableI;qdaTableI;sourceTableI},{[];qdaNoteI;sourceNoteI})

%% Figure 2K - Average Velocity for All 4 Labels

data.exploitPosterior(:) = NaN;
data.sensePosterior(:) = NaN;
data.onPatch(:) = false;
data.exclude(:) = false;

for i = 1:numWorms
    indData = find(data.wormNum == wormNums(i),1);
    indEncount = find(encounter.wormNum == wormNums(i) & ~isnan(encounter.timeEnter) & ...
        strcmp(encounter.expName,expName));
    for j = 1:length(indEncount)
        indTime = encounter.enter(indEncount(j)):encounter.exit(indEncount(j));
        data.exploitPosterior(indData(1) + indTime - 1) = encounter.exploitPosterior(indEncount(j));
        data.sensePosterior(indData(1) + indTime - 1) = encounter.sensePosterior(indEncount(j));
        data.onPatch(indData(1) + indTime - 1) = ...
            strcmp(encounter.label(indEncount(j)),'exploit') | ...
            strcmp(encounter.label(indEncount(j)),'sample') | ...
            strcmp(encounter.label(indEncount(j)),'searchOn');
    end
end

for j = 1:height(info)
    data.exclude(data.plateNum == info.plateNum(j)) = info.exclude(j);
end

indData = (ismember(data.wormNum,wormNums).*~data.exclude.*~data.nearArena) == 1;
weights = indData.*[~data.onPatch,...
    data.onPatch.*(1 - data.exploitPosterior).*(1 - data.sensePosterior),...
    data.onPatch.*(1 - data.exploitPosterior).*data.sensePosterior,...
    data.onPatch.*data.exploitPosterior];
weights(isnan(weights)) = 0;

velocityWeightedMean = sum(data.velocitySmooth.*weights,1,'omitnan')./sum(weights,1,'omitnan')
velocityWeightedStd = sqrt(sum(((data.velocitySmooth - velocityWeightedMean).^2).*weights,1,'omitnan')./...
    ((sum(indData)-1).*sum(weights,1,'omitnan')));

labelNames = {'searchOff','searchOn','sample','exploit'};
velocityBins = 0:450;
kDense = nan(length(labelNames),length(velocityBins));
for i = 1:length(labelNames)
    kDense(i,:) = kde(data.velocitySmooth,'EvaluationPoints',velocityBins,'Bandwidth',400/20,...
        'Weight',weights(:,i));
end

figure; hold on;
for i = 1:length(labelNames)
    plot(velocityBins,kDense(i,:),'LineWidth',1,'Color',eval([labelNames{i},'Color']));
end
set(gca,'XLim',[0 400],'YTick',[],'Units','inches','FontName','Arial',...
    'FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 1 0.5],'Ylim',[0 1.1*max(kDense,[],'all')]);
xlabel('velocity (μm/s)','FontSize',8);
exportgraphics(gca,[saveDir,'fig2k.pdf'])

% Source Data
summaryTableK = table(velocityWeightedMean',velocityWeightedStd','VariableNames',...
    {'weighted mean','weighted std'},'RowNames',{'search (off)','search (on)','sample (on)','exploit (on)'});
summaryNoteK = {'* mean and std of velocity (μm/s) across all time points weighted by probabilities of classification as:';
    '   search (off): 1 - p(on patch)';
    '   search (on): p(on patch)*p(explore)*p(non-sense)';...
    '   sample (on): p(on patch)*p(explore)*p(sense)';...
    '   exploit (on): p(on patch)*p(exploit)*p(sense)'};
kdeTableK = array2table([velocityBins',kDense'],'VariableNames',...
    {'velocity (μm/s)','search (off)','search (on)','sample (on)','exploit (on)'});
kdeNoteK = {'* weighted kernel density estimates (KDE) of distributions of velocity (μm/s)'};
writeSourceData([saveDir,'fig2_data9.xlsx'],{'Figure 2 - Source Data 9';'Figure 2K - Velocity during search, sample, and exploit behaviors'},...
    {'Summary Statistics';'Kernel Density Estimate'},{summaryTableK;kdeTableK},{summaryNoteK;kdeNoteK})

%% Figure 2L - Encounters labeled for all worms

encounterColors = encounter.exploitPosterior.*exploitColor + ...
            (1-encounter.exploitPosterior).*encounter.sensePosterior.*sampleColor + ...
            (1-encounter.exploitPosterior).*(1-encounter.sensePosterior).*searchOnColor;
wormsPerGroup = histcounts(wormGroupInd,(0:numGroups)+0.5)
figure('Position',[500 500 560 800]);
for k = 1:numGroups
    subplot(numGroups,1,1)
    w = wormNums(wormGroupInd == k);
    for i = 1:length(w)
        ind = find(encounter.wormNum == w(i) & strcmp(encounter.expName,expName) & ...
            ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
            strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn')));
        for j = 1:length(ind)
            patch([encounter.timeEnter([ind(j),ind(j)]);encounter.timeExit([ind(j),ind(j)])]./60,...
                i + [-1 0 0 -1]',encounterColors(ind(j),:),'EdgeAlpha',0);
        end
    end
    set(gca,'Units','inches','Visible','on','XTick',[],...
        'FontName','Arial','FontSize',8,'Box','on',...
        'GridLineWidth',0.75,'XLim',[0 60],'YLim',[0 length(w)],...
        'YTick',length(w)/2,'YTickLabel',length(w),'TickDir','out')
    subplotHeight = 0.35*length(w)/max(wormsPerGroup);
    ax = gca; set(gca,'Position',[ax.Position(1), ...
        (numGroups-k+1)*4.5/numGroups + (0.35-subplotHeight)/2, 2, subplotHeight])
end
xlabel('time (min)','FontSize',8); ylabel('# worms','FontSize',8)
set(gca,'XTick',0:10:60)
exportgraphics(gcf,[saveDir,'fig2l.pdf'],'ContentType','vector')

%% Source data
summaryTableL = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',[ones(10,1);12;48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',arrayfun(@(g) sum(indEncounter & G == conditionG(g)),1:numGroups)',...
    'VariableNames',{'label','OD600','growth time (hr)','# worms','# encounters'});
sourceTableL = table(conditionLabels(conditionList),wormList,encounterList,...
    encounter.timeEnter(indEncounter)./60,...
    encounter.duration(indEncounter)./60,log10(encounter.duration(indEncounter)./60),...
    encounter.velocityOn(indEncounter),log10(encounter.velocityOn(indEncounter)),...
    encounter.velocityOnMin(indEncounter),encounter.decelerate(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    encounter.exploitPosterior(indEncounter),1-encounter.exploitPosterior(indEncounter),...
    encounter.sensePosterior(indEncounter),1-encounter.sensePosterior(indEncounter),...
    'VariableNames',{'label','worm #','encounter #','encounter start (min)','encounter duration (min)',...
    'log10(duration)','mean velocity on patch (μm/s)','log10(velocity)',...
    'min. velocity on patch (μm/s)','Δ velocity (μm/s^2)','max. Δ velocity (μm/s^2)',...
    'p(exploit)','p(explore)','p(sense)','p(non-sense)'});
sourceNoteL = {'* posterior probability of classification as explore or exploit as estimated by the Gaussian Mixture Model described in methods and Figure 2 - supplement 6';...
    '* posterior probability of classification as sense or non-sense as estimated by semi-supervised Quadratic Discriminant Analysis (QDA) described in methods, Figure 2I, Figure 2 - supplement 7, and Video 5';...
    '* RGB value calculation: p(exploit)*[0 114 178] + p(explore)*p(sense)*[0 158 115] + p(explore)*p(non-sense)*[230 159 0]'};
writeSourceData([saveDir,'fig2_data10.xlsx'],{'Figure 2 - Source Data 10';'Figure 2L - Encounter classification as explore or exploit and sense or non-sense'},...
    {'Summary Statistics';'Source Data'},{summaryTableL;sourceTableL},{[];sourceNoteL})

%% Figure 2M - p(encounter type) as a function of time

timeBins = ((0:3601)-0.5)./60;
pState = nan(length(timeBins)-1,4,numGroups);
for i = 1:numGroups
    indData = ismember(data.wormNum,wormNums(wormGroupInd == i));
    [~,~,binID] = histcounts(data.timeOffset(indData)./60,timeBins);
    pState(:,:,i) = splitapply(@(p) mean(p,1,'omitnan'),weights(indData,:),binID);
end
pEncounterType = pState(:,2:4,:)./sum(pState(:,2:4,:),2);

smoothSize = 300;
pEncounterTypeSmooth = movmean(pEncounterType(:,:,:),smoothSize,1,'omitnan');
labelNames = {'searchOn','sample','exploit'};
oneHour = (0:3600)./60;

figure('Position',[500 500 560 800]);
for i = 1:numGroups
    subplot(numGroups,1,1); hold on
    for j = 1:3
        patch([oneHour,flip(oneHour)],[sum(pEncounterTypeSmooth(:,1:j,i),2);flip(sum(pEncounterTypeSmooth(:,1:j-1,i),2))]',...
            eval([labelNames{j},'Color']),'EdgeAlpha',0)
    end
    set(gca,'Units','inches','Visible','on','XTick',[],'YTick',[],...
        'FontName','Arial','FontSize',8,'Box','on','TickDir','out',...
        'GridLineWidth',0.75,'XLim',[0 60],'YLim',[0 1],'YAxisLocation','right')
    ax = gca; set(gca,'Position',[ax.Position(1), (12-i+1)*4.5/numGroups, 0.975, 0.35])
    if i == 1
        yticks([0 1]); ylabel('p(encounter)','FontSize',8); set(get(gca,'YLabel'),'Rotation',-90)
    end
end
xlabel('time (min)','FontSize',8); set(gca,'XTick',0:20:60)
exportgraphics(gcf,[saveDir,'fig2m.pdf'],'ContentType','vector')

% Source data
summaryTableM = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',[ones(10,1);12;48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',reshape(mean(pEncounterType(:,1,:),1,'omitnan'),numGroups,1),...
    reshape(mean(pEncounterType(:,2,:),1,'omitnan'),numGroups,1),reshape(mean(pEncounterType(:,3,:),1,'omitnan'),numGroups,1),...
    reshape(mean(sum(pState(:,2:4,:),2),1),numGroups,1),...
    'VariableNames',{'label','OD600','growth time (hr)','# worms',...
    'p(search|on patch)','p(sample|on patch)','p(exploit|on patch)','p(on patch)'})
summaryNoteM = {'* mean probabilities of having an encounter of each type (search, sample, and exploit) are given';...
    '* mean probability of being on patch is also shown'};
sourceTableM = array2table([reshape(repmat(conditionLabels,1,length(oneHour))',[],1),...
    repmat(oneHour',12,1),reshape(permute(pEncounterType,[2,1,3]),3,[])',...
    reshape(sum(pState(:,2:4,:),2),[],1)],...
    'VariableNames',{'label','time (min)','p(search|on patch)','p(sample|on patch)','p(exploit|on patch)','p(on patch)'});
sourceNoteM = {'* probabilities of having an encounter of each type (search, sample, and exploit) are given';...
    '* probability of being on patch is also shown';...
    '* probability values shown in Figure 2M are smoothed using a moving mean filter for improved legibility'};
writeSourceData([saveDir,'fig2_data11.xlsx'],{'Figure 2 - Source Data 11';'Figure 2M - Probability of an encounter type over time'},...
    {'Summary Statistics';'Source Data'},{summaryTableM;sourceTableM},{summaryNoteM;sourceNoteM})

%% Figure 2N,O - # Encounters and Time before an exploit

numEncounters = nan(numWorms,1);
timeExploit = nan(numWorms,1);
doesExploit = nan(numWorms,1);
% onlySearch = nan(numWorms,1);
for i = 1:numWorms
    ind = indEncounter & encounter.wormNum == wormNums(i);
    indExploit = find(ind & strcmp(encounter.label,'exploit'),1,'first');
    indSample = find(ind & ~strcmp(encounter.label,'exploit'));
    doesExploit(i) = ~isempty(indExploit);
    % onlySearch(i) = sum(ind & (strcmp(encounter.label,'exploit') | strcmp(encounter.label,'sample'))) == 0 & ...
    %     sum(ind & strcmp(encounter.label,'searchOn')) > 0;
    if ~doesExploit(i)
        numEncounters(i) = length(indSample);
        timeExploit(i) = max(data.timeOffset(data.wormNum == wormNums(i)));
    else
        numEncounters(i) = sum(indSample < indExploit);
        timeExploit(i) = encounter.timeEnter(indExploit);
    end
end
% numEncountersMean = splitapply(@(X) mean(X,'omitnan'),numEncounters,wormGroup);
%%
timeExploitMean = splitapply(@(X) mean(X,'omitnan'),timeExploit,wormGroupInd);

figure; hold on
% colorInd = wormGroupInd+(numGroups.*doesExploit);
% colorArray = [min(exploreColor + (1 - exploreColor).*(colorValue-0.1),1);...
%     min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1)];
% colorArray = colorArray(unique(colorInd),:);
% gscatter(wormBorder,timeExploit./60,colorInd,colorArray,'.',5)
gscatter(wormBorder,timeExploit./60,doesExploit,[exploreColor;exploitColor],'.',5)
[f,g] = fit(log10(wormBorder)-log10(borderAmp0)+1,timeExploit./60,'logistic4','Lower',[0 -Inf -Inf 0]);
x = logspace(log10(borderAmp0)-1,log10(max(relativeBorder))+1,100);
plot(x,f(log10(x) - log10(borderAmp0) + 1),'k','LineWidth',1)
scatter(relativeBorder,timeExploitMean./60,20,colorValue.*[1 1 1],'filled','MarkerEdgeColor','k')
xscale('log'); xlim([0.005 1000]); ylim([0 60]); yticks(0:20:60)
xticks([0.01 0.1 1 10 100]); xticklabels([0 0.1 1 10 100]); legend('off')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlabel('relative density','FontSize',8); ylabel({'time before';'exploit (min)'},'FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 0.975 0.75])
exportgraphics(gca,[saveDir,'fig2n.pdf'])

% Kendall's correlation (tau is + for monotonically increasing y)
[tauTimeIncrease,pTimeIncrease] = corr(wormBorder,timeExploit./60,'type','Kendall')

% Source data
summaryTableN = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',...
    [repmat(1,10,1);12;48],relativeBorder,...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,...
    timeExploitMean./60,splitapply(@(X) mean(X,'omitnan'),~doesExploit,wormGroupInd),...
    [{tauTimeIncrease};repmat({[]},numGroups-1,1)],[{pTimeIncrease};repmat({[]},numGroups-1,1)],...
    'VariableNames',{'label','OD600','growth time (hr)','relative density',...
    '# worms','saturation','time before first exploit (min)','p(only exploit)','τ','p value'})
summaryNoteN = {'* mean of mean time before first exploitation (min) for each animal';...
    '* in cases where no exploitation is observed (i.e. p(exploit) < 0.5 for all encounters), the last observed time point is used';...
    '* p(only exploit) gives the fraction of animals that were not observed to exploit during the time window';...
    '* p value calculated using Kendall''s τ rank correlation coefficient (test of monotonicity)'};
sigmoidTableN = table(f.a,f.b,f.c,f.d,g.rsquare,'VariableNames',...
    {'a','b','c','d','r-squared'});
sourceTableN = table((1:numWorms)',conditionLabels(wormGroupInd),timeExploit./60,~doesExploit,...
    'VariableNames',{'worm #','label','time before first exploit (min)','only explores'})
sourceNoteN = {'* time (min) before first exploitation (or last observed time point) for each worm';...
    '* animals that only explored (i.e. p(exploit) < 0.5 for all encounters)'};
writeSourceData([saveDir,'fig2_data12.xlsx'],{'Figure 2 - Source Data 12';'Figure 2N - Time before first exploitation'},...
    {'Summary Statistics';'Sigmoid (y = d + (a-d)/(1 + (x/c)^b)';'Source Data'},{summaryTableN;sigmoidTableN;sourceTableN},{summaryNoteN;[];sourceNoteN})

%%
numEncountersMean = splitapply(@(X) mean(X,'omitnan'),numEncounters,wormGroupInd)

figure; hold on
% colorInd = wormGroupInd+(numGroups.*doesExploit);
% colorArray = [min(exploreColor + (1 - exploreColor).*(colorValue-0.1),1);...
%     min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1)];
% colorArray = colorArray(unique(colorInd),:);
% gscatter(wormBorder,numEncounters,colorInd,colorArray,'.',5)
gscatter(wormBorder,numEncounters,doesExploit,[exploreColor;exploitColor],'.',5)
[f,g] = fit(log10(wormBorder)-log10(borderAmp0)+1,numEncounters,'logistic4','Lower',[0 -Inf -Inf 0]);
x = logspace(log10(borderAmp0)-1,log10(max(relativeBorder))+1,100);
plot(x,f(log10(x) - log10(borderAmp0) + 1),'k','LineWidth',1)
scatter(relativeBorder,numEncountersMean,20,colorValue.*[1 1 1],'filled','MarkerEdgeColor','k')
xscale('log'); xlim([0.005 1000]); ylim([0 50]); yticks(0:15:45)
xticks([0.01 0.1 1 10 100]); xticklabels([0 0.1 1 10 100]); legend('off')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlabel('relative density','FontSize',8); ylabel({'# encounters';'before exploit'},'FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 0.975 0.75])
exportgraphics(gca,[saveDir,'fig2o.pdf'])

% Kendall's correlation (tau is + for monotonically increasing y)
[tauTimeIncrease,pTimeIncrease] = corr(wormBorder,numEncounters,'type','Kendall')

% Source data
summaryTableO = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1','1'}',...
    [repmat(1,10,1);12;48],relativeBorder,...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,...
    numEncountersMean,splitapply(@(X) mean(X,'omitnan'),~doesExploit,wormGroupInd),...
    [{tauTimeIncrease};repmat({[]},numGroups-1,1)],[{pTimeIncrease};repmat({[]},numGroups-1,1)],...
    'VariableNames',{'label','OD600','growth time (hr)','relative density',...
    '# worms','saturation','# encounters before first exploit','p(only exploit)','τ','p value'});
summaryNoteO = {'* mean of # encounters before first exploitation (min) for each animal';...
    '* in cases where no exploitation is observed (i.e. p(exploit) < 0.5 for all encounters), the total # of observed encounters is used';...
    '* p(only exploit) gives the fraction of animals that were not observed to exploit during the time window';...
    '* p value calculated using Kendall''s τ rank correlation coefficient (test of monotonicity)'};
sigmoidTableO = table(f.a,f.b,f.c,f.d,g.rsquare,'VariableNames',...
    {'a','b','c','d','r-squared'});
sourceTableO = table((1:numWorms)',conditionLabels(wormGroupInd),numEncounters,~doesExploit,...
    'VariableNames',{'worm #','label','# encounters before first exploit (min)','only explores'})
sourceNoteO = {'* # encounters before first exploitation (or total observed encounters) for each worm';...
    '* animals that only explored (i.e. p(exploit) < 0.5 for all encounters)'};
writeSourceData([saveDir,'fig2_data13.xlsx'],{'Figure 2 - Source Data 13';'Figure 2O - Number of encounters before first exploitation'},...
    {'Summary Statistics';'Sigmoid (y = d + (a-d)/(1 + (x/c)^b)';'Source Data'},{summaryTableO;sigmoidTableO;sourceTableO},{summaryNoteO;[];sourceNoteO})