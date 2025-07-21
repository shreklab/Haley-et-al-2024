%% Load data into workspace

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'geometricData','geometricData.mat'));
allVars = vertcat(modelVars{:});
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
load(fullfile(path,expName,'experimentInfo.mat'),'info');
saveDir = [path,'figures\Fig4\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
searchOffColor = [240 228 66]./255; % yellow
permuteColor = [204 121 167]./255; % pink
exploreColor = [213 94 0]./255; % red
rhoKColor = [61 107 102]./255; % turquoise
tauSColor = [178 73 36]./255; % rust
rhoHColor = [107 61 90]./255; % marroon
rhoEColor = [69 61 109]./255; % indigo
betaColors = [0.5 0.5 0.5; rhoKColor; tauSColor; rhoHColor; rhoEColor];

%% Get conditions to analyze

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));

conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 0.5 & ...
    ~(strcmp(GID.OD600Label,'0.10') & GID.growthCondition == 48));
numGroups = length(conditionG);

%% Get average estimated amplitude of each condition and normalize to OD = 10 (0.5 uL)

borderAmp = splitapply(@(X) mean(X,'omitnan'),encounter.borderAmplitude,G);
relativeBorder = borderAmp(conditionG);
borderAmp10 = borderAmp(strcmp(GID.expName,'foragingConcentration') & ...
    strcmp(GID.OD600Label,'10.00') & GID.lawnVolume == 0.5);
borderAmp0 = 1e-2; % assign 0 to 0.01
relativeBorder = 10.*relativeBorder./borderAmp10;
relativeBorder(strcmp(GID.OD600Label(conditionG),'0.00')) = borderAmp0;

% Assign each border value to a color
colorValue = min(1 - (log(relativeBorder)*0.09 + 0.4),0.8);
conditionLabels = round([0;relativeBorder(2:end)],1,'significant');

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
% indEncounter = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
%     ~encounter.exclude;
indEncounter = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn'));

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%% Get model vars used for this model

[nReps,nBoot] = size(beta{nModels},[2,3]);
allVars = vertcat(modelVars{:});
allVars = allVars(~allVars.exclude & ismember(allVars.conditionNum,conditionG) & ...
    ismember(allVars.wormNum,wormNums),...
    [{'repNum','conditionNum','wormNum','wormID','encounterID','encounterNum'},...
    predictors,response]);

% bootVars = [];
% for b = 1:nBoot
%     ind = arrayfun(@(id) find(allVars.wormID == id),...
%         bootWorm(:,b),'UniformOutput',false);
%     ind = vertcat(ind{:});
%     bootVars = [bootVars;allVars(ind,:)];
% end

%% Figure 4A - Schematic of GLM

% made in illustrator

%% Figure 4B - Schematic of Probabilistic Inclusion of Encounters

rng(1)
skewness = 4;
nEncounters = 10;
sensePosterior = rand(nEncounters,1).^((((skewness^2)-1)*binornd(1,0.3,[nEncounters,1])+1)/skewness);
sensePosterior = flip(sensePosterior); sensePosterior([2,4,5]) = [0.5;0.8;0.9];
exploitPosterior = round(sensePosterior).*rand(nEncounters,1);
exploitPosterior(2:3) = [0.1;0.2];
sense = binornd(1,repmat(sensePosterior,1,4));
exploit = binornd(1,repmat(exploitPosterior,1,4));

figure; hold on;
for i = 1:nEncounters
    text(-0.5,0.5+i,num2str(i),...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',8)
    patch([0 1 1 0],[0 0 1 1]+i,...
        sensePosterior(i).*sampleColor + (1-sensePosterior(i)).*searchOnColor)
    patch([0 1 1 0]+1.5,[0 0 1 1]+i,...
        exploitPosterior(i).*exploitColor + (1-exploitPosterior(i)).*sampleColor)
    for j = 1:4
        if sense(i,j) == 1
            if j == 4
                patch([0 1 1 0]+1.5*j+8,[0 0 1 1]+i,...
                    exploit(i,j).*exploitColor + (1-exploit(i,j)).*sampleColor)
            else
                patch([0 1 1 0]+1.5*j+7,[0 0 1 1]+i,...
                    exploit(i,j).*exploitColor + (1-exploit(i,j)).*sampleColor)
            end
        end

    end
end
for j = 1:3
    text(1.5*j+7.5,nEncounters+1.5,num2str(j),...
            'VerticalAlignment','middle','HorizontalAlignment','center',...
            'FontName','Arial','FontSize',8)
end
text(1.5*4+7.25,nEncounters+1.5,'...',...
            'VerticalAlignment','middle','HorizontalAlignment','center',...
            'FontName','Arial','FontSize',8)
text(1.5*4+8.5,nEncounters+1.5,num2str(100),...
            'VerticalAlignment','middle','HorizontalAlignment','center',...
            'FontName','Arial','FontSize',8)
set(gca,'YDir','reverse','XLim',[-1 25],'YLim',[-1 25],'Units','inches',...
    'FontName','Arial','FontSize',8,'Box','on','XTick',[],'YTick',[])
ax = gca; set(gca,'Position',[0.1 0.1 4.5 4.5])
exportgraphics(gcf,[saveDir,'fig4b.pdf'],'ContentType','vector')

%%
numEncounterSensed = arrayfun(@(n) sum(ismember(modelVars{n}.conditionNum,conditionG) & ismember(modelVars{n}.wormNum,wormNums)),1:length(modelVars));
[mean(numEncounterSensed),std(numEncounterSensed),prctile(numEncounterSensed,[0 50 100])]

%% Figure 4C - Plot beta values for the comprehensive model

nModels = length(beta);
meanBeta = mean(beta{nModels},[2,3]);
pBeta = nModels*2*mean(beta{nModels}.*sign(meanBeta) < 0,[2,3])

figure; hold on
yline(0,'k')
for m = 1:nModels

    v = Violin({beta{nModels}(m,:,:)},m,'ViolinColor',{[0 0 0]},...
        'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
    v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
    v.EdgeColor = betaColors(m,:); v.WhiskerPlot.Color = betaColors(m,:);
    v.MedianPlot.MarkerEdgeColor = betaColors(m,:); v.MeanPlot.Color = betaColors(m,:);
    v.MeanPlot.LineWidth = 0.75;

    if pBeta(m) < 0.001
        sig = '***';
    elseif pBeta(m) < 0.01
        sig = '**';
    elseif pBeta(m) < 0.05
        sig = '*';
    else
        sig = 'n.s.';
    end
    text(m,4.5,sig,'Color',betaColors(m,:),'VerticalAlignment','bottom',...
        'HorizontalAlignment','center','FontSize',8,'FontName','Arial')
end
ylim([-2.25 5.25]); xlim(0.5 + [0 nModels]); xticks(1:nModels); xticklabels({'1','ρ_k','τ_s','ρ_h','ρ_e'})
set(gca,'Units','inches','FontName','Arial','FontSize',8,...
    'Box','on','TickDir','out');
ylabel('β','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
set(gca,'Position',[ax.Position(1:2) 0.975 1.75])
exportgraphics(gcf,[saveDir,'fig4c.pdf'],'ContentType','vector')

% Source data
summaryTableC = table(mean(beta{bestModel},[2,3]),std(beta{bestModel},[],[2,3]),...
    median(beta{bestModel},[2,3]),prctile(beta{bestModel},2.5,[2,3]),...
    prctile(beta{bestModel},97.5,[2,3]),pBeta,...
    'VariableNames',{'mean','std','median','2.5%','97.5%','p value'},'RowNames',{'β_0','β_k','β_s','β_h','β_e'},...
    'DimensionNames',{'coefficient','Variables'});
summaryNoteC = {'* coefficient values fit by logistic regression across 50,000 replicates as described in methods and Figure 4 - supplement 1';...
    '* p value calculated using two-tailed, one-sample bootstrap hypothesis tests with Bonferonni Correction for multiple comparisons'};
[n1,n2] = ndgrid(1:size(beta{nModels},2),1:size(beta{nModels},3));
sourceTableC = array2table([reshape(n1,[],1),reshape(n2,[],1),reshape(beta{nModels},nModels,[])'],...
    'VariableNames',{'encounter replicate #','bootstrap replicate #','β_0','β_k','β_s','β_h','β_e'})
sourceNoteC = {'* coefficient values for all 50,000 combinations of:';....
    '    100 replicates of probabilistic inclusion of sensed encounters (i.e. included encounters if sense = 1 ~ Bern(p(sense = 1)) where p(sense) was previously estimate in Figure 2I and Figure 2 - supplement 7)';...
    '    500 replicates of hierarchically bootstrapped worms (i.e. resampled all worms with replacement and then included all sensed encounters for those animals)'};
writeSourceData([saveDir,'fig4_data1.xlsx'],{'Figure 4 - Source Data 1';'Figure 4C - Coefficients of linear regression model'},...
    {'Summary Statistics';'Source Data'},{summaryTableC;sourceTableC},{summaryNoteC;sourceNoteC})

%% Figure 4D - Plot p(first exploit) vs. encounter #

condition01 = conditionG(GID.growthCondition(conditionG) == 0 & strcmp(GID.OD600Label(conditionG),'1.00'));
condition10 = conditionG(strcmp(GID.OD600Label(conditionG),'10.00'));

histVars = allVars(ismember(allVars.conditionNum,[condition01,condition10]) & ...
    ismember(allVars.wormNum,wormNums),:);
[wormID,~,wormRep] = unique(histVars(:,{'conditionNum','wormNum','repNum'}));
betaNames = {'β_0','+ β_k ρ_k','+ β_s τ_s','+ β_h ρ_h','+ β_e ρ_e'};
maxEncounter = max(histVars.encounterNum);

% Functions to get histogram of outcomes
histBins = @(c) unique(histVars.encounterNum(histVars.conditionNum == c));
histOutcomes = @(pE,c) splitapply(@mean,pE(histVars.conditionNum == c),histVars.encounterNum(histVars.conditionNum == c));

% Simulate a set of outcomes given the observed probabilities of acceptance
outcomes_observed = binornd(1,repmat(histVars.exploitK,1,1));

% Get encounter # of first accept
firstVisit_observed = nan(height(wormID),1);
[c,ia] = unique(wormRep.*outcomes_observed(:),'stable');
firstVisit_observed(c(c > 0)) = histVars.encounterNum(ia(c > 0));

% Plot observed histograms
figure('Position',[400 400 800 400]);
for m = 1:nModels
    
    % Estimate probability of accept from model (outcome = 1)
    meanBeta = mean(beta{m},[2,3]);
    theseVars = [ones(height(histVars),1),table2array(histVars(:,predictors(1:m-1)))];
    histVars.(['exploitK_GLM',num2str(m)]) = glmval(meanBeta,theseVars,'logit','constant','off');

    % Simulate a set of outcomes given the modeled probabilities of acceptance
    outcomes_simulated = binornd(1,repmat(histVars.(['exploitK_GLM',num2str(m)]),1,1));

    % Plot p(first exploit)
    firstVisit_simulated = nan(height(wormID),1);
    [c,ia] = unique(wormRep.*outcomes_simulated(:),'stable');
    firstVisit_simulated(c(c > 0)) = histVars.encounterNum(ia(c > 0));

    % Plot simulated exploitations
    subplot(2,nModels+1,m); hold on
    h.(['predicted',num2str(m),'_01']) = histogram(firstVisit_simulated(wormID.conditionNum == condition01,:),0.5+(0:maxEncounter),...
        'Normalization','probability','FaceColor',colorValue(condition01).*[1 1 1],'FaceAlpha',1);
    
    xlim([0 maxEncounter]); ylim([0 0.35])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 2.9 0.75 0.75])
    title(betaNames{m},'FontSize',9,'FontWeight','normal','FontName','Times New Roman','FontAngle','italic');
 
    subplot(2,nModels+1,m+nModels+1); hold on
    h.(['predicted',num2str(m),'_10']) = histogram(firstVisit_simulated(wormID.conditionNum == condition10,:),0.5+(0:maxEncounter),...
        'Normalization','probability','FaceColor',colorValue(condition10).*[1 1 1],'FaceAlpha',1);
    xlim([0 maxEncounter]); ylim([0 0.35])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 2 0.75 0.75])
    if m == 1
        xticks(0:10:30); yticks(0:0.1:0.3)
        ylabel('p(first exploit)','FontSize',8); xlabel('encounter #','FontSize',8)
    end

end

% Plot observed first exploit
subplot(2,nModels+1,nModels+1); hold on
h.observed_01 = histogram(firstVisit_observed(wormID.conditionNum == condition01,:),0.5+(0:maxEncounter),...
    'Normalization','probability','FaceColor',colorValue(condition01).*[1 1 1],'FaceAlpha',1);
xlim([0 maxEncounter]); ylim([0 0.35])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
set(gca,'Position',[5.5 2.9 0.75 0.75])

subplot(2,nModels+1,2*(nModels+1)); hold on
h.observed_10 = histogram(firstVisit_observed(wormID.conditionNum == condition10,:),0.5+(0:maxEncounter),...
    'Normalization','probability','FaceColor',colorValue(condition10).*[1 1 1],'FaceAlpha',1);
xlim([0 maxEncounter]); ylim([0 0.35])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
set(gca,'Position',[5.5 2 0.75 0.75])

exportgraphics(gcf,[saveDir,'fig4d.pdf'],'ContentType','vector')

% Source data
meanBeta = cell(nModels,nModels);
for m = 1:nModels
    meanBeta(1:m,m) = num2cell(mean(beta{m},[2,3]));
end
betaTableD = cell2table(meanBeta,'RowNames',{'β_0','β_k','β_s','β_h','β_e'},...
    'VariableNames',{'β_0','+ β_k ρ_k','+ β_s τ_s','+ β_h ρ_h','+ β_e ρ_e'},...
    'DimensionNames',{'coefficient','Variables'});
betaNoteD = {'* coefficient values fit by logistic regression across 50,000 replicates as described in methods and Figure 4 - supplement 1'};

histogramTableD = table([ones(maxEncounter,1);10*ones(maxEncounter,1)],repmat((1:maxEncounter)',2,1),...
    [h.predicted1_01.Values,h.predicted1_10.Values]',[h.predicted2_01.Values,h.predicted2_10.Values]',...
    [h.predicted3_01.Values,h.predicted3_10.Values]',[h.predicted4_01.Values,h.predicted4_10.Values]',...
    [h.predicted5_01.Values,h.predicted5_10.Values]',[h.observed_01.Values,h.observed_10.Values]',...
    'VariableNames',{'label','encounter #','β_0','+ β_k ρ_k','+ β_s τ_s','+ β_h ρ_h','+ β_e ρ_e','observed'});
histogramNoteD = {'* probability of exploitation occuring for the first time as a function of the number of encounters';...
    '* exploitation events were simulated for 100 replicates using the p(exploit) and p(sense) as described in methods, Figure 2H, and Figure 2 - supplement 6'};

[wormID,~] = unique(allVars(:,{'conditionNum','wormNum','encounterNum','repNum'}),'stable');
[~,wormList] = ismember(allVars.wormNum,wormNums);
[~,conditionList] = ismember(allVars.conditionNum,conditionG);
for m = 1:nModels
    meanBeta = mean(beta{m},[2,3]);
    theseVars = [ones(height(allVars),1),table2array(allVars(:,predictors(1:m-1)))];
    allVars.(['exploitK_GLM',num2str(m)]) = glmval(meanBeta,theseVars,'logit','constant','off');
end
sourceTableD = table(conditionLabels(conditionList),wormList,allVars.repNum,...
    allVars.encounterNum,allVars.exploitK_GLM1,allVars.exploitK_GLM2,...
    allVars.exploitK_GLM3,allVars.exploitK_GLM4,allVars.exploitK_GLM5,allVars.exploitK,...
    'VariableNames',{'label','worm #','replicate #','sensed encounter #',...
    'β_0','+ β_k ρ_k','+ β_s τ_s','+ β_h ρ_h','+ β_e ρ_e','observed'});
sourceNoteD = {'* predicted probability of exploitation estimated for each sensed encounter using models with increasing number of covariates';...
    '* observed probability of exploitation estimated using GMM classifier as described in methods, Figure 2H, and Figure 2 - supplement 6'
    '* exploitations were simulated from a Bernoulli distribution of these probabilities'};

writeSourceData([saveDir,'fig4_data2.xlsx'],{'Figure 4 - Source Data 2';...
    'Figure 4D - Probability of exploitation occuring for the first time as a function of the number of encounters'},...
    {'Coefficient values for nested GLMs';'Probability distributions';'Source Data'},...
    {betaTableD;histogramTableD;sourceTableD},{betaNoteD;histogramNoteD;sourceNoteD})