%% Load data into workspace

expName = 'foragingMutants';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'geometricData','geometricData.mat'));
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);

saveDir = [path,'figures\FigS17\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
searchOffColor = [240 228 66]./255; % yellow
permuteColor = [204 121 167]./255; % pink
exploreColor = [213 94 0]./255; % red
rhoKColor = [61 107 102]./255; % turquoise
rhoHColor = [107 61 90]./255; % marroon
rhoEColor = [69 61 109]./255; % indigo
betaColors = [0.5 0.5 0.5; rhoKColor; rhoHColor; rhoEColor];

%% Get conditions to analyze

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));

conditionG = find(strcmp(GID.expName,expName) & strcmp(GID.strainID,'N2'));
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
    ~encounter.exclude;

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%% Figure S17A - Example traces of well-fed and food-deprived

wormNum = [30, 35]; % [FD, WF]
for i = 1:length(wormNum)
    % plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
end

%%
numEncounterTotal = arrayfun(@(g) sum(indEncounter & G == g),conditionG)
numEncounterSensed = arrayfun(@(g) arrayfun(@(n) sum(ismember(modelVars{n}.conditionNum,conditionG(g)) & ...
    ismember(modelVars{n}.wormNum,wormNums)),1:length(modelVars)),1:length(conditionG),'UniformOutput',false);
numEncounterSensed = vertcat(numEncounterSensed{:})
[numEncounterTotal,mean(numEncounterSensed,2),std(numEncounterSensed,[],2),prctile(numEncounterSensed,[0 50 100],2)]

%% Figure S17B - Encounters labeled for all worms

mutantColors = [permuteColor;0 0 0];

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
        'YTick',[0 length(w)],'TickDir','out')
    ylabel('# worms','FontSize',8)
    subplotHeight = 0.85*length(w)/max(wormsPerGroup);
    ax = gca; set(gca,'Position',[ax.Position(1),3-1.1*k, 2, subplotHeight])
    title(GID.strainName{conditionG(k)},'FontWeight','normal','FontSize',8,'FontName','Arial','Color',mutantColors(k,:))
end
xlabel('time (min)','FontSize',8)
set(gca,'XTick',0:10:60)
exportgraphics(gcf,[saveDir,'figS17b.pdf'],'ContentType','vector')

%% Estimate betas without time

nModels = length(beta);
beta_noSatiety = nan(nModels-1,length(modelVars),size(bootWorm,2));
for n = 1:length(modelVars)
    include = ~modelVars{n}.exclude & ismember(modelVars{n}.conditionNum,...
        find(strcmp(GID.expName,'foragingConcentration') & GID.lawnVolume == 0.5 & ...
        ~(strcmp(GID.OD600Label,'0.10') & GID.growthCondition == 48)));
    theseVars = modelVars{n}(include,['wormID',predictors,'exploitK']);
    for b = 1:size(bootWorm,2)
        ind = arrayfun(@(id) find(theseVars.wormID == id),...
            bootWorm(:,b),'UniformOutput',false);
        ind = vertcat(ind{:});
        mdl = fitglm(theseVars(ind,:),'linear','Distribution','binomial',...
            'PredictorVars',predictors(~strcmp(predictors,'timeOffSinceExploit')),...
            'ResponseVar','exploitK','binomialSize',1,...
            'Options',statset('MaxIter',10000),'LikelihoodPenalty','none');
        beta_noSatiety(:,n,b) = mdl.Coefficients.Estimate;
    end
    n
end

%% Figure S17C - Plot beta values without satiety

meanBeta = mean(beta_noSatiety,[2,3]);
pBeta = size(beta_noSatiety,1)*2*mean(beta_noSatiety.*sign(meanBeta) < 0,[2,3])
indBeta = find(~strcmp(predictors,'timeOffSinceExploit'));

figure; hold on
yline(0,'k')
for m = 1:size(beta_noSatiety,1)

    v = Violin({beta_noSatiety(m,:,:)},m,'ViolinColor',{[0 0 0]},...
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
ylim([-2.25 5.25]); xlim(0.5 + [0 nModels-1]); xticks(1:nModels-1); xticklabels({'1','ρ_k','ρ_h','ρ_e'})
set(gca,'Units','inches','FontName','Arial','FontSize',8,...
    'Box','on','TickDir','out');
ylabel('β','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
set(gca,'Position',[ax.Position(1:2) 0.975 1.75])
exportgraphics(gcf,[saveDir,'figS17c.pdf'],'ContentType','vector')

%% Figure S17D - Plot p(first exploit) for WF and FD (predicted + observed)

allVars = vertcat(modelVars{:});
allVars = allVars(ismember(allVars.conditionNum,conditionG) & ismember(allVars.wormNum,wormNums),:);
[wormID,~,wormRep] = unique(allVars(:,{'conditionNum','wormNum','repNum'}));
betaNames = {{'estimated','without satiety'},{'estimated','with satiety'}};
betaNames2 = {'without satiety','with satiety'};
maxEncounter = max(allVars.encounterNum);
nModels = length(beta);
nReps = 1;
rng(1)
mutantColors = [permuteColor;colorValue(2)*ones(1,3)];

% Add starved time to timeOffSinceExploit
starvedDuration = arrayfun(@(w) mean(hours(info.starvedDuration(cellfun(@(worms) any(ismember(worms,w)),info.wormNum)))),allVars.wormNum);
allVars.timeOffSinceExploit = allVars.timeOffSinceExploit + starvedDuration;

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
figure('Position',[400 0 1200 1200]);
for m = 1:2
    
    % Estimate probability of accept from model (outcome = 1)        
    switch m
        case 2
            meanBeta = mean(beta{end},[2,3]);
            theseVars = [ones(height(allVars),1),table2array(allVars(:,predictors))];
            allVars.(['exploitK_GLM',num2str(m)]) = glmval(meanBeta,theseVars,'logit','constant','off');
        case 1
            meanBeta = mean(beta_noSatiety,[2,3]);
            theseVars = [ones(height(allVars),1),table2array(allVars(:,predictors(~strcmp(predictors,'timeOffSinceExploit'))))];
            allVars.(['exploitK_GLM',num2str(m)]) = glmval(meanBeta,theseVars,'logit','constant','off');
    end

    % Simulate a set of outcomes given the modeled probabilities of acceptance
    outcomes_simulated = binornd(1,repmat(allVars.(['exploitK_GLM',num2str(m)]),1,nReps));

    % Plot p(first exploit)
    firstVisit_simulated = nan(height(wormID),nReps);
    for i = 1:nReps
        [c,ia] = unique(wormRep.*outcomes_simulated(:,i),'stable');
        firstVisit_simulated(c(c > 0),i) = allVars.encounterNum(ia(c > 0));
    end
    for g = 1:numGroups
        thisCondition = conditionG(g);

        % Plot simulated exploitations
        subplot(numGroups,nModels+1,numGroups*(nModels+1)); hold on
        h = histogram(firstVisit_simulated(wormID.conditionNum == thisCondition,:),0.5+(0:maxEncounter),...
            'Normalization','probability','FaceColor',mutantColors(g,:),'FaceAlpha',1);

        xlim([0 maxEncounter]); ylim([0 1])
        set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
        set(gca,'Position',[m*0.9 (numGroups-g+1)*0.8 0.75 0.75])
        if g == 1
            title(betaNames{m},'FontSize',8,'FontWeight','normal','FontName','Arial');
        end

        if g == numGroups && m == 1
            xticks(0:3:9); yticks(0:0.5:1)
            ylabel('p(first exploit)','FontSize',8); xlabel('encounter #','FontSize',8)
        end
    end
end

% Plot observed first exploit
for g = 1:numGroups
    thisCondition = conditionG(g);
    subplot(2,nModels+1,nModels+1); hold on
    h = histogram(firstVisit_observed(wormID.conditionNum == thisCondition,:),0.5+(0:maxEncounter),...
        'Normalization','probability','FaceColor',mutantColors(g,:),'FaceAlpha',1);
    xlim([0 maxEncounter]); ylim([0 1])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[2.8 (numGroups-g+1)*0.8 0.75 0.75])
    if g == 1
        title('observed','FontSize',8,'FontWeight','normal','FontName','Arial');
    end
end
exportgraphics(gcf,[saveDir,'figS17d.pdf'],'ContentType','vector')