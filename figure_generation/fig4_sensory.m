%% Load data into workspace

expName = 'foragingSensory';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'geometricData','geometricData.mat'),'modelVars','predictors');
model.WF = load(fullfile(path,'geometricData','geometricData_well-fed.mat'));
model.OSM6 = load(fullfile(path,'geometricData','geometricData_osm-6.mat'));
model.MEC4 = load(fullfile(path,'geometricData','geometricData_mec-4.mat'));
lambda.WF = load(fullfile(path,'geometricData','geometricData_well-fed_ridge.mat'));
lambda.OSM6 = load(fullfile(path,'geometricData','geometricData_osm-6_ridge.mat'));
lambda.MEC4 = load(fullfile(path,'geometricData','geometricData_mec-4_ridge.mat'));
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,'midpoint.mat'),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
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
strainColor = [0 0 0; permuteColor; exploreColor];

%% Get conditions to analyze

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));

strainNames = {'well-fed','mec-4','osm-6'};
strainID = {'WF','MEC4','OSM6'};
conditionG = find(ismember(GID.strainName,strainNames) & ...
        (strcmp(GID.expName,expName) | strcmp(GID.expName,'foragingMutants')));
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

%% Get ids for worms used in this model

[wormID,wormIDInfo] = findgroups(table(G,encounter.wormNum,'VariableNames',{'conditionNum','wormNum'}));

% Remove worms tracked for less than 75% of the video
wormNumsSensory = unique(wormID(ismember(G,conditionG) & ~encounter.exclude & ...
    ismember(G,find(strcmp(GID.expName,expName)))));
framesTrackedSensory = arrayfun(@(w) sum(~data.noTrack(data.wormNum == w))/...
    sum(info.numFrames(info.plateNum == unique(data.plateNum(data.wormNum == w)))), wormIDInfo.wormNum(wormNumsSensory));
wormNumsMutants = unique(wormID(ismember(G,conditionG) & ~encounter.exclude & ...
    ismember(G,find(strcmp(GID.expName,'foragingMutants')))));
load(fullfile(path,'foragingMutants','experimentInfo.mat'),'info');
load(fullfile(path,'foragingMutants','midpoint.mat'),'data');
framesTrackedMutants = arrayfun(@(w) sum(~data.noTrack(data.wormNum == w))/...
    sum(info.numFrames(info.plateNum == unique(data.plateNum(data.wormNum == w)))), wormIDInfo.wormNum(wormNumsMutants));
wormNums = [wormNumsMutants;wormNumsSensory];
framesTracked = [framesTrackedMutants;framesTrackedSensory];
wormNums = wormNums(framesTracked >= 0.75);
numWorms = length(wormNums);

% Get condition id of each worm
wormGroup = arrayfun(@(w) unique(G(wormID == w & ismember(G,conditionG))),wormNums);
[~,~,wormGroupInd] = unique(wormGroup); % continuous group id

% Get indices of encounters
indEncounter = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & ~strcmp(encounter.label,'') & ~strcmp(encounter.label,'searchOff');

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(wormID == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%%
numEncounterTotal = arrayfun(@(g) sum(indEncounter & G == g),conditionG)
numEncounterSensed = arrayfun(@(n) sum(ismember(modelVars{n}.conditionNum,conditionG) & ismember(modelVars{n}.wormID,wormNums)),1:length(modelVars));
[sum(numEncounterTotal),mean(numEncounterSensed),std(numEncounterSensed),prctile(numEncounterSensed,[0 50 100])]

%% Figure 4G - Plot p(first exploit) vs. encounter #

allVars = vertcat(modelVars{:});
allVars = allVars(ismember(allVars.conditionNum,conditionG) & ...
    ismember(allVars.wormID,wormNums) & ~allVars.exclude,:);
[wormID,~,wormRep] = unique(allVars(:,{'conditionNum','wormNum','repNum'}));
maxEncounter = max(allVars.encounterNum);
densities = {'1.00','10.00'};
nReps = 1;

% Functions to get histogram of outcomes
histBins = @(c) unique(allVars.encounterNum(allVars.conditionNum == c));
histOutcomes = @(pE,c) splitapply(@mean,pE(allVars.conditionNum == c),allVars.encounterNum(allVars.conditionNum == c));

% Simulate a set of outcomes given the observed probabilities of acceptance
outcomes_observed = binornd(1,repmat(allVars.exploitK,1,nReps));

% Get encounter # of first accept
firstVisit_observed = nan(height(wormID),nReps);
for i = 1:nReps
    [c,ia] = unique(wormRep.*outcomes_observed(:,i),'stable');
    firstVisit_observed(c(c > 0),i) = allVars.encounterNum(ia(c > 0));
end

% Plot observed histograms
figure('Position',[400 400 800 400]);
for s = 1:length(strainNames)

    for d = 1:length(densities)
    conditionID = find(strcmp(GID.strainName(conditionG),strainNames{s}) & ...
        strcmp(GID.OD600Label(conditionG),densities{d}));

    % Plot observed first exploit
    subplot(length(densities),length(strainNames),length(densities)*length(strainNames)); hold on
    plotColor = mean(min(strainColor(s,:) + (1 - strainColor(s,:)).*(colorValue(conditionID)-0.1),1),1);
    h.([strainID{s},'_',num2str(d)]) = histogram(firstVisit_observed(ismember(wormID.conditionNum,conditionG(conditionID)),:),0.5+(0:maxEncounter),...
        'Normalization','probability','FaceColor',plotColor,'FaceAlpha',1);
    xlim([0 maxEncounter]); % ylim([0 0.35])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on',...
        'TickDir','out','XTick',[],'YLim',[0 0.9],'YTick',[]);
    set(gca,'Position',[0.9*s 4-0.9*d 0.75 0.75])
    if d == 1
        title(strainNames{s},'FontSize',8,'FontWeight','normal','FontName','Arial');
    end
    if d == length(densities) && s == 1
        xticks(0:5:15); yticks(0:0.3:0.9)
        ylabel('p(first exploit)','FontSize',8)
        xlabel('encounter #','FontSize',8)
    end
    end
end
exportgraphics(gcf,[saveDir,'fig4g.pdf'],'ContentType','vector')

%% Source data

histogramTableD = table([ones(maxEncounter,1);10*ones(maxEncounter,1)],repmat((1:maxEncounter)',2,1),...
    [h.WF_1.Values,h.WF_2.Values]',[h.MEC4_1.Values,h.MEC4_2.Values]',[h.OSM6_1.Values,h.OSM6_2.Values]',...
    'VariableNames',{'label','encounter #','well-fed','mec-4','osm-6'});
histogramNoteD = {'* probability of exploitation occuring for the first time as a function of the number of encounters';...
    '* exploitation events were simulated for 100 replicates using the p(exploit) and p(sense) as described in methods, Figure 2H, and Figure 2 - supplement 6'};

[wormID,~] = unique(allVars(:,{'conditionNum','wormNum','encounterNum','repNum'}),'stable');
[~,wormList] = ismember(allVars.wormID,wormNums);
[~,conditionList] = ismember(allVars.conditionNum,conditionG);
sourceTableD = table(GID.strainName(conditionG(conditionList)),allVars.ODK,...
    wormList,allVars.repNum,allVars.encounterNum,allVars.exploitK,...
    'VariableNames',{'strain','label','worm #','replicate #','sensed encounter #','p(exploit)'});
sourceNoteD = {'* observed probability of exploitation estimated using GMM classifier as described in methods, Figure 2H, and Figure 2 - supplement 6'
    '* exploitations were simulated from a Bernoulli distribution of these probabilities'};
sourceTableD = sortrows(sourceTableD,{'replicate #','worm #'});

writeSourceData([saveDir,'fig4_data4.xlsx'],{'Figure 4 - Source Data 4';...
    'Figure 4G - Probability of exploitation occuring for the first time as a function of the number of encounters for animals with sensory mutations'},...
    {'Probability distributions';'Source Data'},...
    {histogramTableD;sourceTableD},{histogramNoteD;sourceNoteD})

%% Figure 4H - Plot beta values for each strain

% Compute p-value as mean of differences b/w mutants and wild-type
% distributions of covariates
nModels = length(predictors) + 1;
for s = 1:length(strainNames)
    betaAll(s,:,:) = reshape(model.(strainID{s}).beta,nModels,[]);
end
meanBeta = mean(betaAll,3)
mu = mean(betaAll,3);
sigma = std(betaAll,[],3);
Z = (mu(2:end,:) - mu(1,:))./(sigma(2:end,:).^2 + sigma(1,:).^2).^(0.5);
pTwoTail = 2*min(normcdf(Z),normcdf(Z,'upper')); % two-tailed
pOneTail = (sign(mu(1,:)) < 0) + sign(mu(1,:)).*normcdf(Z); % one-tailed (in the direction of 0 from WT)
p = [pTwoTail(:,1),pOneTail(:,2),pTwoTail(:,3),pOneTail(:,4:5)]; % no hypothesis about beta_0 or tau_s, hypothesize rhos are less salient
[~,pAdj] = benjaminiHochberg(reshape(p',1,[]),0.05);
pAdj = reshape(pAdj',nModels,[])';

for s = 1:length(strainNames)
    beta = model.(strainID{s}).beta;
    figure; hold on
    yline(0,'k')
    for m = 1:nModels
        v = Violin({beta(m,:,:)},m,'ViolinColor',{strainColor(s,:)},...
            'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
        v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
        v.EdgeColor = betaColors(m,:); v.WhiskerPlot.Color = betaColors(m,:);
        v.MedianPlot.MarkerEdgeColor = betaColors(m,:); v.MeanPlot.Color = betaColors(m,:);
        v.MeanPlot.LineWidth = 0.75;

        if s > 1
            pBeta = pAdj(s-1,:);
            if pBeta(m) < 0.001
                sig = '***';
            elseif pBeta(m) < 0.01
                sig = '**';
            elseif pBeta(m) < 0.05
                sig = '*';
            elseif p(s-1,m) < 0.05
                sig = '†';
            else
                sig = 'n.s.';
            end
        text(m,7,sig,'Color',betaColors(m,:),'VerticalAlignment','bottom',...
            'HorizontalAlignment','center','FontSize',8,'FontName','Arial')
        end
    end
    ylim([-3.25 8]); xlim([0.5 5.5]); xticks(1:5); xticklabels({'1','ρ_k','τ_s','ρ_h','ρ_e'})
    set(gca,'Units','inches','FontName','Arial','FontSize',8,...
        'Box','on','TickDir','out');
    ylabel('β','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
    ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
    set(gca,'Position',[ax.Position(1:2) 0.975 1.75])
    exportgraphics(gcf,[saveDir,'fig4h_',strainNames{s},'.pdf'],'ContentType','vector')
end

%% Source data

modelLambda = nan(size(strainNames));
for s = 1:length(strainNames)
    meanLogLike = mean(lambda.(strainID{s}).logLike,[2,3]);
    plot(lambda.(strainID{s}).lambda,meanLogLike,'k')
    [maxLogLike,bestLambda] = max(meanLogLike);
    modelLambda(s) = lambda.(strainID{s}).lambda(bestLambda);
end

summaryTableH_WF = table(mean(betaAll(1,:,:),3)',std(betaAll(1,:,:),[],3)',...
    median(betaAll(1,:,:),3)',prctile(betaAll(1,:,:),2.5,3)',...
    prctile(betaAll(1,:,:),97.5,3)',[modelLambda(1);repmat({[]},length(predictors),1)],...
    'VariableNames',{'mean','std','median','2.5%','97.5%','λ'},'RowNames',{'β_0','β_k','β_s','β_h','β_e'},...
    'DimensionNames',{'coefficient','Variables'});
summaryNoteH_WF = {'* coefficient values fit by logistic regression across 50,000 replicates as described in methods, Figure 4C, and Figure 4 - supplement 1';...
    '* regression models were regularized using a ridge parameter λ which was optimized to maximize the mean log-likelihood'};

summaryTableH_MEC4 = table(mean(betaAll(2,:,:),3)',std(betaAll(2,:,:),[],3)',...
    median(betaAll(2,:,:),3)',prctile(betaAll(2,:,:),2.5,3)',...
    prctile(betaAll(2,:,:),97.5,3)',[modelLambda(2);repmat({[]},length(predictors),1)],...
    p(1,:)',pAdj(1,:)',...
    'VariableNames',{'mean','std','median','2.5%','97.5%','λ','p-value','adjusted p-value'},'RowNames',{'β_0','β_k','β_s','β_h','β_e'},...
    'DimensionNames',{'coefficient','Variables'});
summaryNoteH_MEC4 = {'* coefficient values fit by logistic regression across 50,000 replicates as described in methods, Figure 4C, and Figure 4 - supplement 1';...
 '* regression models were regularized using a ridge parameter λ which was optimized to maximize the mean log-likelihood';...
 '* p value calculated using two-sample mean of differences between mec-4 and N2 distribution with Benjamini-Hochberg correction for multiple comparisons';...
 '* one-tailed test for density coefficents β_k, β_h, and β_e (alternative hypothesis = less influence of density in sensory mutants)';...
 '* two-tailed test for other coefficients β_0 and β_s (alternative hypothesis = altered influence in sensory mutants)'};

summaryTableH_OSM6 = table(mean(betaAll(3,:,:),3)',std(betaAll(3,:,:),[],3)',...
    median(betaAll(3,:,:),3)',prctile(betaAll(3,:,:),2.5,3)',...
    prctile(betaAll(3,:,:),97.5,3)',[modelLambda(3);repmat({[]},length(predictors),1)],...
    p(2,:)',pAdj(2,:)',...
    'VariableNames',{'mean','std','median','2.5%','97.5%','λ','p-value','adjusted p-value'},'RowNames',{'β_0','β_k','β_s','β_h','β_e'},...
    'DimensionNames',{'coefficient','Variables'});
summaryNoteH_OSM6 = {'* coefficient values fit by logistic regression across 50,000 replicates as described in methods, Figure 4C, and Figure 4 - supplement 1';...
 '* regression models were regularized using a ridge parameter λ which was optimized to maximize the mean log-likelihood';...
 '* p value calculated using two-sample mean of differences between osm-6 and N2 distribution with Benjamini-Hochberg correction for multiple comparisons';...
 '* one-tailed test for density coefficents β_k, β_h, and β_e (alternative hypothesis = less influence of density in sensory mutants)';...
 '* two-tailed test for other coefficients β_0 and β_s (alternative hypothesis = altered influence in sensory mutants)'};

[n1,n2] = ndgrid(1:size(beta,2),1:size(beta,3));
sourceTableH_WF = array2table([reshape(n1,[],1),reshape(n2,[],1),permute(betaAll(1,:,:),[2,3,1])'],...
    'VariableNames',{'encounter replicate #','bootstrap replicate #','β_0','β_k','β_s','β_h','β_e'});
sourceTableH_MEC4 = array2table([reshape(n1,[],1),reshape(n2,[],1),permute(betaAll(3,:,:),[2,3,1])'],...
    'VariableNames',{'encounter replicate #','bootstrap replicate #','β_0','β_k','β_s','β_h','β_e'});
sourceTableH_OSM6 = array2table([reshape(n1,[],1),reshape(n2,[],1),permute(betaAll(3,:,:),[2,3,1])'],...
    'VariableNames',{'encounter replicate #','bootstrap replicate #','β_0','β_k','β_s','β_h','β_e'});
sourceTableH = addvars([sourceTableH_WF;sourceTableH_MEC4;sourceTableH_OSM6],...
    reshape(repmat(strainNames,size(betaAll,3),1),[],1),...
    'NewVariableNames','strain','Before','encounter replicate #')
sourceNoteH = {'* coefficient values for all 50,000 combinations of:';....
    '    100 replicates of probabilistic inclusion of sensed encounters (i.e. included encounters if sense = 1 ~ Bern(p(sense = 1)) where p(sense) was previously estimate in Figure 2I and Figure 2 - supplement 7)';...
    '    500 replicates of hierarchically bootstrapped worms (i.e. resampled all worms with replacement and then included all sensed encounters for those animals)'};

writeSourceData([saveDir,'fig4_data5.xlsx'],{'Figure 4 - Source Data 5';'Figure 4H - Coefficients of linear regression model'},...
    {'Summary Statistics (N2)';'Summary Statistics (mec-4)';'Summary Statistics (osm-6)';'Source Data'},...
    {summaryTableH_WF;summaryTableH_MEC4;summaryTableH_OSM6;sourceTableH},...
    {summaryNoteH_WF;summaryNoteH_MEC4;summaryNoteH_OSM6;sourceNoteH})