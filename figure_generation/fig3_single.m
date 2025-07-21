%% Load data into workspace

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
thresh = readtable([path,'foragingMini\encounterThresholds.csv']);
saveDir = [path,'figures\Fig3\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
searchOffColor = [240 228 66]./255; % yellow
exploreColor = [213 94 0]./255; % red

%% Get condition ids

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume','growthCondition',...
    'lawnOD600','peptone'}));
conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 20);
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

mean(vertcat(info.lawnDiameterMean{ismember(info.plateNum,unique(encounter.plateNum(indEncounter)))}))

%% Figure 3A - Concentrations tested w/ effective OD600 (Single)

figure; hold on;
xValues = 1:numGroups; xValues(end) = xValues(end) + 0.5;
bar(xValues,relativeBorder,'FaceColor','flat','CData',colorValue.*[1 1 1])
for i = 1:numGroups
    t = text(xValues(i),relativeBorder(i),num2str(round(relativeBorder(i),1,'significant')),...
        'FontName','Arial','FontSize',6,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
    if i == 1
        t.String = '0'; t.Position = [1 0.05 0];
    end
end
text(0.75,300,'Ø 8.3 mm','HorizontalAlignment','left','VerticalAlignment','middle',...
    'FontSize',8,'FontName','Arial')
xlim([0.25 numGroups+1.25])
xticks(xValues); xticklabels({'0','0.05','0.1','0.5','1','2','3','4','5','10','1'})
yscale('log'); ylim([0.05 1600]); yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 2*(numGroups+1)/13 0.75])
ylabel('relative density','FontSize',8)
exportgraphics(gca,[saveDir,'fig3a.pdf'])

% Source data
relativeBorder_std = splitapply(@(X) std(X,'omitnan'),10*encounter.borderAmplitude./borderAmp10,G);
relativeBorder_std = relativeBorder_std(conditionG);
summaryTableA = table({'0','0.05','0.1','0.5','1','2','3','4','5','10','1'}',...
    [ones(10,1);48],relativeBorder,relativeBorder_std,conditionLabels,colorValue,...
    'VariableNames',{'OD600','growth time (hr)','mean relative density','std relative density','label','saturation'})
summaryNoteA = {'* relative density was estimated as described in methods and Figure 2 - supplement 1-3';...
    '* bacteria-free patches were assigned a relative density 0.01 for use in log-based calculations';...
    '* condition labels use rounded relative density estimates for ease of reading';...
    '* saturation gives the grayscale color value from 0 (black) to 1 (white) for each condition'};
writeSourceData([saveDir,'fig3_data1.xlsx'],{'Figure 3 - Source Data 1';'Figure 3A - Relative density of large (20 μL) bacterial patches'},...
    {'Summary Statistics'},{summaryTableA;},{summaryNoteA})

%% Figure 3D - Example Traces (Single)

wormNum = 514;
%plotTracks(info,data,wormNum,'timeOffset',saveDir,bodyPart,'yes');

wormNum = [501,502];
for i = 1:length(wormNum)
    plotTracks(info,data,wormNum(i),'timeOffset',[path,'figures\Fig4\'],bodyPart,'yes');
end

%% Figure 3F - Encounters labeled for all worms (single)

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
    subplotHeight = 0.41*length(w)/40;
    ax = gca; set(gca,'Position',[ax.Position(1), ...
        (numGroups-k+1)*2.5/numGroups + (0.41-subplotHeight)/2, 2, subplotHeight])
end
xlabel('time (min)','FontSize',8); ylabel('# worms','FontSize',8)
set(gca,'XTick',0:10:60)
exportgraphics(gcf,[saveDir,'fig3f.pdf'],'ContentType','vector')

% Source data
[~,wormList] = ismember(encounter.wormNum(indEncounter),wormNums);
[~,conditionList] = ismember(G(indEncounter),conditionG);
encounterList = arrayfun(@(w) 1:sum(encounter.wormNum(indEncounter)==w),unique(encounter.wormNum(indEncounter)),'UniformOutput',false);
encounterList = [encounterList{:}]';
summaryTableF = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1'}',[ones(10,1);48],...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',arrayfun(@(g) sum(indEncounter & G == conditionG(g)),1:numGroups)',...
    'VariableNames',{'label','OD600','growth time (hr)','# worms','# encounters'});
sourceTableF = table(conditionLabels(conditionList),wormList,encounterList,...
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
sourceNoteF = {'* posterior probability of classification as explore or exploit as estimated by the Gaussian Mixture Model described in methods, Figure 2H, and Figure 2 - supplement 6';...
    '* posterior probability of classification as sense or non-sense as estimated by semi-supervised Quadratic Discriminant Analysis (QDA) described in methods, Figure 2I, Figure 2 - supplement 7, and Video 5';...
    '* RGB value calculation: p(exploit)*[0 114 178] + p(explore)*p(sense)*[0 158 115] + p(explore)*p(non-sense)*[230 159 0]'};
writeSourceData([saveDir,'fig3_data3.xlsx'],{'Figure 3 - Source Data 3';'Figure 3F - Encounter classification as explore or exploit and sense or non-sense'},...
    {'Summary Statistics';'Source Data'},{summaryTableF;sourceTableF},{[];sourceNoteF})

%% Figure 3H - Time before an exploit (single)

numEncounters = nan(numWorms,1);
timeExploit = nan(numWorms,1);
doesExploit = nan(numWorms,1);
for i = 1:numWorms
    ind = indEncounter & encounter.wormNum == wormNums(i);
    indExploit = find(ind & strcmp(encounter.label,'exploit'),1,'first');
    indSample = find(ind & ~strcmp(encounter.label,'exploit'));
    doesExploit(i) = ~isempty(indExploit);
    if ~doesExploit(i)
        numEncounters(i) = length(indSample);
        timeExploit(i) = max(data.timeOffset(data.wormNum == wormNums(i)));
    else
        numEncounters(i) = sum(indSample < indExploit);
        timeExploit(i) = encounter.timeEnter(indExploit);
    end
end
timeExploitMean = splitapply(@(X) mean(X,'omitnan'),timeExploit,wormGroupInd);

figure; hold on
% colorInd = wormGroupInd+(numGroups.*doesExploit);
% colorArray = [max(min(exploreColor + (1 - exploreColor).*(colorValue-0.1),1),0);...
%     max(min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1),0)];
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
xlabel('relative density','FontSize',8); ylabel({'time before','exploit (min)'},'FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 0.975 0.75])
exportgraphics(gca,[saveDir,'fig3h.pdf'])

% Kendall's correlation (tau is + for monotonically increasing y)
[tauTimeIncrease,pTimeIncrease] = corr(wormBorder,timeExploit,'type','Kendall')

% Source data
summaryTableH = table(conditionLabels,...
    {'0','0.05','0.1','0.5','1','2','3','4','5','10','1'}',...
    [ones(10,1);48],relativeBorder,...
    arrayfun(@(g) sum(wormGroupInd == g),1:numGroups)',colorValue,...
    timeExploitMean./60,splitapply(@(X) mean(X,'omitnan'),~doesExploit,wormGroupInd),...
    [{tauTimeIncrease};repmat({[]},numGroups-1,1)],[{pTimeIncrease};repmat({[]},numGroups-1,1)],...
    'VariableNames',{'label','OD600','growth time (hr)','relative density',...
    '# worms','saturation','time before first exploit (min)','p(only exploit)','τ','p value'})
summaryNoteH = {'* mean of mean time before first exploitation (min) for each animal';...
    '* in cases where no exploitation is observed (i.e. p(exploit) < 0.5 for all encounters), the last observed time point is used';...
    '* p(only exploit) gives the fraction of animals that were not observed to exploit during the time window';...
    '* p value calculated using Kendall''s τ rank correlation coefficient (test of monotonicity)'};
sigmoidTableH = table(f.a,f.b,f.c,f.d,g.rsquare,'VariableNames',...
    {'a','b','c','d','r-squared'});
sourceTableH = table((1:numWorms)',conditionLabels(wormGroupInd),timeExploit./60,~doesExploit,...
    'VariableNames',{'worm #','label','time before first exploit (min)','only explores'})
sourceNoteH = {'* time (min) before first exploitation (or last observed time point) for each worm';...
    '* animals that only explored (i.e. p(exploit) < 0.5 for all encounters)'};
writeSourceData([saveDir,'fig3_data5.xlsx'],{'Figure 3 - Source Data 5';'Figure 3H - Time before first exploitation'},...
    {'Summary Statistics';'Sigmoid (y = d + (a-d)/(1 + (x/c)^b)';'Source Data'},{summaryTableH;sigmoidTableH;sourceTableH},{summaryNoteH;[];sourceNoteH})
%%

figure; hold on
colorInd = wormGroupInd+(numGroups.*doesExploit);
colorArray = [max(min(exploreColor + (1 - exploreColor).*(colorValue-0.1),1),0);...
    max(min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1),0)];
colorArray = colorArray(unique(colorInd),:);
gscatter(wormBorder,numEncounters,colorInd,colorArray,'.',5)
% f = fit(log10(wormBorder)-log10(borderAmp0)+1,numEncounters,'logistic4','Lower',[0 -Inf -Inf 0],'Upper',[15 0 Inf 15]);
% x = logspace(log10(borderAmp0)-1,log10(max(relativeBorder))+1,100);
% plot(x,f(log10(x) - log10(borderAmp0) + 1),'k','LineWidth',1)
xscale('log'); xlim([0.005 1000]); ylim([0 15]); yticks(0:5:15)
xticks([0.01 0.1 1 10 100]); xticklabels([0 0.1 1 10 100]); legend('off')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlabel('relative density','FontSize',8); ylabel({'# encounters';'before exploit'},'FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 0.975 0.75])
exportgraphics(gca,[saveDir,'fig3h.pdf'])

% Kendall's correlation (tau is + for monotonically increasing y)
[tauTimeIncrease,pTimeIncrease] = corr(wormBorder,numEncounters,'type','Kendall')