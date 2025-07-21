%% Load data into workspace

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'geometricData','geometricData.mat'));
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
load(fullfile(path,expName,'experimentInfo.mat'),'info');

saveDir = [path,'figures\FigS16\'];

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

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%% Figure S16 - Plot p(first exploit) vs. encounter #

allVars = vertcat(modelVars{:});
allVars = allVars(ismember(allVars.conditionNum,conditionG),:);
[wormID,~,wormRep] = unique(allVars(:,{'conditionNum','wormNum','repNum'}));
betaNames = {'β_0','+ β_k ρ_k','+ β_s τ_s','+ β_h ρ_h','+ β_e ρ_e'};
maxEncounter = max(allVars.encounterNum);
nModels = length(beta);
nReps = 1;
rng(11)

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

[~,cOrder] = sort(colorValue,'descend');

% Plot observed histograms
KLD = nan(nModels,length(unique(allVars.repNum)),numGroups);
figure('Position',[400 0 1200 1200]);
for m = 1:nModels

    % Estimate probability of accept from model (outcome = 1)
    meanBeta = mean(beta{m},[2,3]);
    theseVars = [ones(height(allVars),1),table2array(allVars(:,predictors(1:m-1)))];
    allVars.(['exploitK_GLM',num2str(m)]) = glmval(meanBeta,theseVars,'logit','constant','off');

    % Simulate a set of outcomes given the modeled probabilities of acceptance
    outcomes_simulated = binornd(1,repmat(allVars.(['exploitK_GLM',num2str(m)]),1,nReps));

    % Plot p(first exploit)
    firstVisit_simulated = nan(height(wormID),nReps);
    for i = 1:nReps
        [c,ia] = unique(wormRep.*outcomes_simulated(:,i),'stable');
        firstVisit_simulated(c(c > 0),i) = allVars.encounterNum(ia(c > 0));
    end
    for g = 1:numGroups
        thisCondition = conditionG(cOrder(g));

        % Plot simulated exploitations
        subplot(numGroups,nModels+1,numGroups*(nModels+1)); hold on
        h = histogram(firstVisit_simulated(wormID.conditionNum == thisCondition,:),0.5+(0:maxEncounter),...
            'Normalization','probability','FaceColor',colorValue(g).*[1 1 1],'FaceAlpha',1);

        xlim([0 maxEncounter]); ylim([0 0.9])
        set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
        set(gca,'Position',[m*0.9 (numGroups-g+1)*0.8 0.75 0.75])
        if g == 1
            title(betaNames{m},'FontSize',9,'FontWeight','normal','FontName','Times New Roman','FontAngle','italic');
        end

        if g == numGroups && m == 1
            xticks(0:10:30); yticks(0:0.3:0.9)
            ylabel('p(first exploit)','FontSize',8); xlabel('encounter #','FontSize',8)
        end

        P = arrayfun(@(n) histcounts(firstVisit_observed(wormID.conditionNum == thisCondition & wormID.repNum == n,:),...
            0.5+(0:maxEncounter),'Normalization','probability'),unique(wormID.repNum),'UniformOutput',false);
        Q = arrayfun(@(n) histcounts(firstVisit_simulated(wormID.conditionNum == thisCondition & wormID.repNum == n,:),...
            0.5+(0:maxEncounter),'Normalization','probability'),unique(wormID.repNum),'UniformOutput',false);
        P = vertcat(P{:}); Q = vertcat(Q{:});
        P = max(P,1e-10); Q = max(Q,1e-10);
        KLD(m,:,g) = sum(P.*log(P./Q),2);
    end
end

% Plot observed first exploit
for g = 1:numGroups
    thisCondition = conditionG(cOrder(g));
    subplot('Position',[1 0 0.1 0.1]); hold on
    h = histogram(firstVisit_observed(wormID.conditionNum == thisCondition,:),0.5+(0:maxEncounter),...
        'Normalization','probability','FaceColor',colorValue(g).*[1 1 1],'FaceAlpha',1);
    xlim([0 maxEncounter]); ylim([0 0.9])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[5.5 (numGroups-g+1)*0.8 0.75 0.75])

    subplot('Position',[1 0 0.1 0.1]); hold on
    for m = 1:nModels
        v = Violin({KLD(m,:,g)},m,'ViolinColor',{colorValue(g).*[1 1 1]},...
            'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
        v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
    end
    xlim(0.5+[0 nModels]); ylim([0 max(KLD,[],'all')])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[6.7 (numGroups-g+1)*0.8 2.25 0.75])

    if g == numGroups
        set(gca,'XTick',1:5,'YTick',0:3:9,'XTickLabels',betaNames)
        ylabel('KL Divergence','FontName','Arial','FontSize',8)
        ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 8;
    end
end

exportgraphics(gcf,[saveDir,'figS16.pdf'],'ContentType','vector')