%% Load data into workspace

expName = 'foragingSensory';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
model.WF = load(fullfile(path,'geometricData','geometricData_well-fed.mat'));
model.OSM6 = load(fullfile(path,'geometricData','geometricData_osm-6.mat'));
model.MEC4 = load(fullfile(path,'geometricData','geometricData_mec-4.mat'));
lambda.WF = load(fullfile(path,'geometricData','geometricData_well-fed_ridge.mat'));
lambda.OSM6 = load(fullfile(path,'geometricData','geometricData_osm-6_ridge.mat'));
lambda.MEC4 = load(fullfile(path,'geometricData','geometricData_mec-4_ridge.mat'));
load(fullfile(path,'geometricData','geometricData.mat'),'modelVars','predictors');
load(fullfile(path,'geometricData','geometricData_well-fed_subsample.mat'),'beta');
saveDir = [path,'figures\FigS19\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
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
GID.lawnOD600 = cellfun(@(OD) mean(str2num(OD)),GID.OD600Label);
[G2,GID2] = findgroups(GID(conditionG,{'strainName','lawnOD600'}));
numGroups = height(GID2);
groupOrder = [7:9,1:6];

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
wormGroupInd = arrayfun(@(g) ismember(wormGroup,conditionG(G2 == g)),groupOrder,'UniformOutput',false)';
[~,wormGroupInd] = max([wormGroupInd{:}],[],2); % continuous group id
% [~,~,wormGroupInd] = unique(wormGroup); % continuous group id

% Get indices of encounters
indEncounter = ismember(wormID,wormNums) & ~encounter.exclude & ~strcmp(encounter.label,'') & ~strcmp(encounter.label,'searchOff');

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(wormID == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%%  Figure S19A - Plot example traces

wormNum = [79,89,87,93,98,108,103,73,122]; % ['WF 1', 'WF 5', 'WF 10', 'MEC-4 1', 'MEC-4 5', 'MEC-4 10','OSM-6 1', 'OSM-6 5', 'OSM-6 10']
% for i = 1:length(wormNum)
%     plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
% end

%% Figure S19B - Encounters labeled for all worms

encounterColors = encounter.exploitPosterior.*exploitColor + ...
            (1-encounter.exploitPosterior).*encounter.sensePosterior.*sampleColor + ...
            (1-encounter.exploitPosterior).*(1-encounter.sensePosterior).*searchOnColor;
wormsPerGroup = histcounts(wormGroupInd,(0:numGroups)+0.5)
plotColor = nan(numGroups,3);
figure('Position',[500 500 560 800]);
for k = 1:numGroups
    subplot(numGroups,1,1)
    w = wormNums(wormGroupInd == k);
    for i = 1:length(w)
        ind = find(wormID == w(i) & ...
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
        'YTick',length(w)/2,'YTickLabel',length(w),'TickDir','in')
    subplotHeight = 0.7*length(w)/max(wormsPerGroup);
    ax = gca; set(gca,'Position',[ax.Position(1), ...
        (numGroups-k+1)*5.5/numGroups + (0.7-subplotHeight)/2, 2, subplotHeight])
    s = strcmp(strainNames,GID2.strainName(groupOrder(k)));
    plotColor(k,:) = mean(min(strainColor(s,:) + (1 - strainColor(s,:)).*...
        (colorValue(ismember(conditionG,unique(wormGroup(wormGroupInd == k))))-0.1),1),1);
end
xlabel('time (min)','FontSize',8)
set(gca,'XTick',0:10:60)

subplot('Position',[0.6 0.1 0.1 0.1])
image(permute(reshape(plotColor',[3 3 3]),[2,3,1]))
exportgraphics(gcf,[saveDir,'figS19b.pdf'],'ContentType','vector')

%% Figure S19C - Plot p(first exploit) vs. encounter #

allVars = vertcat(modelVars{:});
allVars = allVars(ismember(allVars.conditionNum,conditionG) & ...
    ismember(allVars.wormID,wormNums) & ~allVars.exclude,:);
[wormID,~,wormRep] = unique(allVars(:,{'conditionNum','wormNum','repNum'}));
maxEncounter = max(allVars.encounterNum);
densities = {'1.00','5.00','10.00'};
strainColor = [0 0 0; permuteColor; exploreColor];
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
    histogram(firstVisit_observed(ismember(wormID.conditionNum,conditionG(conditionID)),:),0.5+(0:maxEncounter),...
        'Normalization','probability','FaceColor',plotColor,'FaceAlpha',1);
    xlim([0 maxEncounter]); % ylim([0 0.35])
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on',...
        'TickDir','out','XTick',[],'YLim',[0 1],'YTick',[]);
    set(gca,'Position',[0.9*s 4-0.9*d 0.75 0.75])
    if d == 1
        title(strainNames{s},'FontSize',8,'FontWeight','normal','FontName','Arial');
    end
    if d == length(densities) && s == 1
        xticks(0:5:15); yticks(0:0.5:1)
        ylabel('p(first exploit)','FontSize',8)
        xlabel('encounter #','FontSize',8)
    end
    end
end
exportgraphics(gcf,[saveDir,'figS19c.pdf'],'ContentType','vector')

%% Figure S19D - Plot null beta values for each strain

nModels = length(predictors)+1;
for s = 1:length(strainNames)
    beta = model.(strainID{s}).betaShuffle;
    figure; hold on
    yline(0,'k')
    for m = 1:nModels
        v = Violin({beta(m,:,:)},m,'ViolinColor',{strainColor(s,:)},...
            'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
        v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
        v.EdgeColor = betaColors(m,:); v.WhiskerPlot.Color = betaColors(m,:);
        v.MedianPlot.MarkerEdgeColor = betaColors(m,:); v.MeanPlot.Color = betaColors(m,:);
        v.MeanPlot.LineWidth = 0.75;
    end
    ylim([-5 5]); xlim([0.5 5.5]); yticks(-4:2:4); xticks(1:5); xticklabels({'1','ρ_k','τ_s','ρ_h','ρ_e'})
    set(gca,'Units','inches','FontName','Arial','FontSize',8,...
        'Box','on','TickDir','out');
    ylabel('β^*','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
    ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
    set(gca,'Position',[ax.Position(1:2) 0.975 1.75])
    exportgraphics(gcf,[saveDir,'figS19d_',strainNames{s},'.pdf'],'ContentType','vector')
end

%% Figure S19E - Ridge regression parameter optimization

figure
for s = 1:length(strainNames)
    subplot(1,length(strainNames),s); hold on
    meanLogLike = mean(lambda.(strainID{s}).logLike,[2,3]);
    plot(lambda.(strainID{s}).lambda,meanLogLike,'k')
    [maxLogLike,bestLambda] = max(meanLogLike);
    modelLambda = lambda.(strainID{s}).lambda(bestLambda);
    scatter(modelLambda,maxLogLike,'k','filled')
    logLikeRange = [min(meanLogLike) max(meanLogLike)];
    text(modelLambda,maxLogLike+diff(logLikeRange)*0.1,['λ = ',num2str(modelLambda)],...
        'VerticalAlignment','bottom','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',8)
    xticks(lambda.(strainID{s}).lambda([1,end]))
    yticks(round(logLikeRange,4))
    xlim(lambda.(strainID{s}).lambda([1,end]))
    ylim(logLikeRange + diff(logLikeRange).*[-0.1 0.5])
    xlabel('λ','FontName','Times New Roman','FontSize',8)
    set(gca,'Units','inches','FontName','Arial','FontSize',8,...
        'Box','on','TickDir','out');
    if s == 1
        ylabel('log-likelihood','FontName','Arial','FontSize',8)
    end
    ax = gca;
    set(gca,'Position',[ax.Position(1:2) 0.975 0.75])
end
exportgraphics(gcf,[saveDir,'figS19e.pdf'],'ContentType','vector')

%% Figure S19D - Subsamples of foragingConcentration have non-significant rho_h

nSample = size(beta,4);
nModels = size(beta,1);
betaNames = {'β_0','β_k','β_s','β_h','β_e'};
rhoKColor = [61 107 102]./255; % turquoise
tauSColor = [178 73 36]./255; % rust
rhoHColor = [107 61 90]./255; % marroon
rhoEColor = [69 61 109]./255; % indigo
betaColors = [0 0 0; rhoKColor; tauSColor; rhoHColor; rhoEColor];
figure
for m = 1:nModels
    subplot('Position',[0.9 0.9 0.1 0.1]); hold on
    yline(0,'k')
    for s = 1:nSample
        p = mean(beta(m,:,:,s) > 0,'all');
        if p < 0.05 || p > 0.95
            vColor = [0 0 0];
        else
            vColor = permuteColor;
        end
        v = Violin({beta(m,:,:,s)},s,'ViolinColor',{vColor},...
            'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
        v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
    end
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on',...
        'TickDir','out','XTick',[],'XLim',[0 nSample+1]);
    ylabel(betaNames{m},'FontSize',9,'FontName','Times New Roman',...
        'FontWeight','bold','FontAngle','italic','Color',betaColors(m,:))
    set(gca,'Position',[0.5 3-0.5*m 4 0.45])
end
set(gca,'XTick',0:10:nSample); xlabel('subsample #','FontSize',8)

exportgraphics(gcf,[saveDir,'figS19d.pdf'],'ContentType','vector')