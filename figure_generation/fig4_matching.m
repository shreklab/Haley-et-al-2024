%% Load data into workspace

expName = 'foragingMatching';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'geometricData','geometricData.mat'));
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
saveDir = [path,'figures\Fig4\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
searchOffColor = [240 228 66]./255; % yellow
permuteColor = [204 121 167]./255; % pink
exploreColor = [213 94 0]./255; % red

%% Get conditions to analyze

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));

conditionG = find(strcmp(GID.expName,expName));
numGroups = length(conditionG);

%% Get average estimated amplitude of each condition and normalize to OD = 10 (0.5 uL)

borderAmp = splitapply(@(X) mean(X,'omitnan'),encounter.borderAmplitude,G);
relativeBorder = borderAmp(conditionG);
borderAmp10 = borderAmp(strcmp(GID.expName,'foragingConcentration') & ...
    strcmp(GID.OD600Label,'10.00') & GID.lawnVolume == 0.5);
borderAmp0 = 1e-2; % assign 0 to 0.01
relativeBorder = 10.*relativeBorder./borderAmp10;
relativeBorder(strcmp(GID.OD600Label(conditionG),'0.00')) = borderAmp0;

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
    ~encounter.exclude & ~strcmp(encounter.label,'') & ~strcmp(encounter.label,'searchOff');

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%%  Figure 4E - Plot example trace

% made in illustrator

%%
numEncounterTotal = arrayfun(@(g) sum(indEncounter & G == g),conditionG)
numEncounterSensed = arrayfun(@(n) sum(ismember(modelVars{n}.conditionNum,conditionG) & ismember(modelVars{n}.wormNum,wormNums)),1:length(modelVars));
[sum(numEncounterTotal),mean(numEncounterSensed),std(numEncounterSensed),prctile(numEncounterSensed,[0 50 100])]

%% Group encounters by ODK + ODH and ODK + ODE

% Concatenate model replicates
allVars = vertcat(modelVars{:});
allVars = allVars(ismember(allVars.conditionNum,conditionG) & ~allVars.exclude & ...
    ismember(allVars.wormNum,wormNums),:);

% Compute log likelihood manually (why does matlab use the gammaln offsets?)
% y = observed p(exploit), n = # binomial trials, p = predicted p(exploit)
% logLikelihood = @(y,n,p) sum(y.*log(p) + (1-y).*log(1-p) + gammaln(n+1) - gammaln(y+1) - gammaln(n-y+1));
logLikelihood = @(y,n,p) sum(y.*log(p) + (1-y).*log(1-p));
logLikeRatio = @(llComplex,llSimple) -2*(llSimple - llComplex); % logLikelihood
logLikeRatioTest = @(LLR) chi2cdf(LLR,1,'upper'); % logLikelihoodRatio

% Calculate estimates of p(exploit) for each nested model
nModels = length(beta);
logLike = nan(nModels,1);
for m = 1:nModels
    meanBeta = mean(beta{m},[2,3]);
    theseVars = [ones(height(allVars),1),table2array(allVars(:,predictors(1:m-1)))];
    allVars.(['exploitK_GLM',num2str(m)]) = glmval(meanBeta,theseVars,'logit','constant','off');
    logLike(m) = logLikelihood(allVars.exploitK,1,allVars.(['exploitK_GLM',num2str(m)]));
end

% Calculate likelihood ratio (w/ and w/o history)
LLR = logLikeRatio(logLike(end),logLike(3))
pLLR = logLikeRatioTest(LLR)

% Get groups of pairings of ODK + ODH and ODK + ODE
[pairID.ODH,pairInfo.ODH] = findgroups(allVars(:,{'ODK','ODH'}));
[pairID.ODE,pairInfo.ODE] = findgroups(allVars(:,{'ODK','ODE'}));

histType = {'ODH','ODE'};

% Get group stats
densityK_tick = arrayfun(@(OD) mean(allVars.densityK(allVars.ODK == OD),'omitnan'),unique(allVars.ODK));
densityK_min = min(allVars.densityK);
densityK_max = max(allVars.densityK);
timeOffSinceExploit_mean = mean(allVars.timeOffSinceExploit);
densityH_mean = mean(allVars.densityH);
densityE_mean = mean(allVars.densityE);
for h = 1:length(histType)
    pairInfo.(histType{h}).densityK_observed = splitapply(@mean,allVars.densityK,pairID.(histType{h}));
    pairInfo.(histType{h}).densityH_observed = splitapply(@mean,allVars.densityH,pairID.(histType{h}));
    pairInfo.(histType{h}).densityE_observed = splitapply(@mean,allVars.densityE,pairID.(histType{h}));
    pairInfo.(histType{h}).pE_observed = splitapply(@mean,allVars.exploitK,pairID.(histType{h}));
    pairInfo.(histType{h}).pE_predicted = splitapply(@mean,allVars.(['exploitK_GLM',num2str(nModels)]),pairID.(histType{h}));
    pairInfo.(histType{h}).pE_predicted_noH = splitapply(@mean,allVars.exploitK_GLM3,pairID.(histType{h}));
    pairInfo.(histType{h}).timeOff = splitapply(@mean,allVars.timeOffSinceExploit,pairID.(histType{h}));
end

%% Figure 4F - Plot p(exploit) for combinations of ODK and ODH

% colors = interp1([0;1],[sampleColor;exploitColor],linspace(0,1,100));
pE_range = [0 1];
colors = jet;
modelType = {'pE_predicted_noH','pE_predicted','pE_observed'};
modelTitles = {{'estimated','no history'},{'estimated','with history'},'observed'};
histLabels = {'ρ_h','ρ_e'};
modelNum = [3,5];

for h = 1:length(histType)
    if strcmp(histType{h},'ODH')
        rho = 'densityH_observed';
    elseif strcmp(histType{h},'ODE')
        rho = 'densityE_observed';
    end
    % pE_range = prctile([pairInfo.(histType{h}).pE_observed;...
    %     pairInfo.(histType{h}).pE_predicted;...
    %     pairInfo.(histType{h}).pE_predicted_noH],[0 100]);

    primaryVertices = 10.^[pairInfo.(histType{h}).densityK_observed(pairInfo.(histType{h}).(histType{h}) < 200),...
        pairInfo.(histType{h}).(rho)(pairInfo.(histType{h}).(histType{h}) < 200)];
    primaryFaces = [1 2 5 4; 4 5 8 7; 2 3 6 5; 5 6 9 8];
    axLimits = prctile(primaryVertices,[0 100],'all');

    extraVertices = 10.^[densityK_min densityK_min;...
        pairInfo.(histType{h}).densityK_observed(pairInfo.(histType{h}).ODK == 5 & pairInfo.(histType{h}).(histType{h}) == 1) densityK_min;...
        densityK_max densityK_min;...
        densityK_min pairInfo.(histType{h}).(rho)(pairInfo.(histType{h}).ODK == 1 & pairInfo.(histType{h}).(histType{h}) == 5);...
        densityK_max pairInfo.(histType{h}).(rho)(pairInfo.(histType{h}).ODK == 10 & pairInfo.(histType{h}).(histType{h}) == 5);...
        densityK_min densityK_max;...
        pairInfo.(histType{h}).densityK_observed(pairInfo.(histType{h}).ODK == 5 & pairInfo.(histType{h}).(histType{h}) == 10) densityK_max;...
        densityK_max densityK_max];
    extraFaces = [10 11 4 1; 11 12 7 4; 12 14 8 7; 14 17 9 8; 17 16 6 9;...
        16 15 3 6; 15 13 2 3; 13 10 1 2];
    
    figure;
    for m = 1:length(modelType)
        
        primaryColors = [pairInfo.(histType{h}).(modelType{m})(pairInfo.(histType{h}).(histType{h}) < 200)];
        extraColors = [pairInfo.(histType{h}).(modelType{m})(pairInfo.(histType{h}).(histType{h}) == 1);...
            pairInfo.(histType{h}).(modelType{m})(pairInfo.(histType{h}).(histType{h}) == 5 & pairInfo.(histType{h}).ODK ~= 5);...
            pairInfo.(histType{h}).(modelType{m})(pairInfo.(histType{h}).(histType{h}) == 10)];
        subplot(1,3,m); hold on

        [X,Y] = ndgrid(linspace(axLimits(1),axLimits(2),100),linspace(axLimits(1),axLimits(2),100));
        switch h
            case 1
                meshVars = [ones(10000,1),log10(X(:)),repmat(timeOffSinceExploit_mean,10000,1),...
                    log10(Y(:)),repmat(densityE_mean,10000,1)];
            case 2
                meshVars = [ones(10000,1),log10(X(:)),repmat(timeOffSinceExploit_mean,10000,1),...
                    repmat(densityH_mean,10000,1),log10(Y(:))];
        end

         switch m
            case 1
                meanBeta = mean(beta{3},[2,3]);
                historyMesh = glmval(meanBeta,meshVars(:,1:length(meanBeta)),'logit','constant','off');
                imagesc('XData',axLimits,'YData',axLimits,'CData',reshape(historyMesh,size(X))');
                switch h
                    case 1
                        scatterVars = [ones(9,1),log10(primaryVertices(:,1)),repmat(timeOffSinceExploit_mean,9,1),...
                            log10(primaryVertices(:,2)),repmat(densityE_mean,9,1)];
                    case 2
                        scatterVars = [ones(9,1),log10(primaryVertices(:,1)),repmat(timeOffSinceExploit_mean,9,1),...
                            repmat(densityH_mean,9,1),log10(primaryVertices(:,2))];
                end
                historyScatter = glmval(meanBeta,scatterVars(:,1:length(meanBeta)),'logit','constant','off');
                scatter(primaryVertices(:,1),primaryVertices(:,2),...
                     'filled','CData',historyScatter,'MarkerEdgeColor','k','SizeData',100)
            case 2
                meanBeta = mean(beta{5},[2,3]);
                historyMesh = glmval(meanBeta,meshVars,'logit','constant','off');
                imagesc('XData',axLimits,'YData',axLimits,'CData',reshape(historyMesh,size(X))');
                switch h
                    case 1
                        scatterVars = [ones(9,1),log10(primaryVertices(:,1)),repmat(timeOffSinceExploit_mean,9,1),...
                            log10(primaryVertices(:,2)),repmat(densityE_mean,9,1)];
                    case 2
                        scatterVars = [ones(9,1),log10(primaryVertices(:,1)),repmat(timeOffSinceExploit_mean,9,1),...
                            repmat(densityH_mean,9,1),log10(primaryVertices(:,2))];
                end
                historyScatter = glmval(meanBeta,scatterVars(:,1:length(meanBeta)),'logit','constant','off');
                scatter(primaryVertices(:,1),primaryVertices(:,2),...
                     'filled','CData',historyScatter,'MarkerEdgeColor','k','SizeData',100)
             case 3
                 patch('Vertices',[primaryVertices;extraVertices],...
                     'Faces',[primaryFaces;extraFaces],...
                     'FaceVertexCData',[primaryColors;extraColors],...
                     'FaceColor','interp','EdgeColor','none')
                 scatter(primaryVertices(:,1),primaryVertices(:,2),...
                     'filled','CData',primaryColors,'MarkerEdgeColor','k','SizeData',100)
        end
       
       
        colormap(colors); clim(pE_range);
        xlim(axLimits); ylim(axLimits)
        set(gca,'YDir','normal','TickDir','out','XTick',10.^densityK_tick,...
            'XTickLabel',[1 4 7],'YTick',10.^densityK_tick,'YTickLabel',[1 4 7],...
            'Units','inches','FontName','Arial','FontSize',8);
        if m < length(modelType)
            set(gca,'Position',[0.5+1.1*(m-1) 2 0.75 0.75])
        else
            set(gca,'Position',[0.5+1.15*(m-1) 2 0.75 0.75])
        end
        if m == 1
            xlabel('10^{ρ_k}','FontSize',9,'FontName','Times New Roman','FontAngle','italic');
            ylabel(['10^{',histLabels{h},'}'],'FontSize',9,'FontName','Times New Roman','FontAngle','italic');
            ax = gca; ax.XRuler.TickLabelGapOffset = 5; ax.YRuler.TickLabelGapOffset = 5;
        else
            set(gca,'XTick',[],'YTick',[])
        end
        t = title(modelTitles{m},'FontSize',8,'FontWeight','normal');
        t.Position(2) = t.Position(2)+0.6;
    end

    colormap(colors); c = colorbar('Units','inches');
    set(c,'Position',[3.8 2 0.1 0.75],'Limits',pE_range,'Ticks',0:1);
    set(c.Label,'String','p(exploit)','FontSize',8,'Rotation',-90)
    exportgraphics(gcf,[saveDir,'fig4f_',histType{h},'.pdf'])
end

%% Source data

meanBeta = cell(nModels,nModels);
for m = 1:nModels
    meanBeta(1:m,m) = num2cell(mean(beta{m},[2,3]));
end
betaTableF = cell2table(meanBeta(:,[3,5]),'RowNames',{'β_0','β_k','β_s','β_h','β_e'},...
    'VariableNames',{'no history','with history'},...
    'DimensionNames',{'coefficient','Variables'});
betaNoteF = {'* coefficient values fit by logistic regression across 50,000 replicates as described in methods and Figure 4 - supplement 1';...
    '* model without history: β_0 + β_k ρ_k + β_s τ_s';...
    '* model with history: β_0 + β_k ρ_k + β_s τ_s + β_h ρ_h + β_e ρ_e'};

logLikeTableF = table(logLike(3),logLike(5),LLR,pLLR,'VariableNames',...
    {'no history','with history','likelihood-ratio test','p'});
logLikeNoteF = {'* log-likelihood = sum(y.*log(p) + (1-y).*log(1-p))'};

summaryTableF_ODH = array2table([ceil(0.7*pairInfo.ODH{:,{'ODK','ODH'}}),...
    pairInfo.ODH{:,{'densityK_observed','densityH_observed','pE_predicted_noH','pE_predicted','pE_observed'}}],...
    'VariableNames',{'label (10^ρ_k)','label (10^ρ_h)','ρ_k','ρ_h','no history','with history','observed'});
summaryTableF_ODH(pairInfo.ODH.ODH == 200,:) = [];
summaryTableF_ODH = addvars(summaryTableF_ODH,repmat(timeOffSinceExploit_mean,9,1),...
    'NewVariableNames',{'τ_s'},'After','ρ_k')
summaryTableF_ODH = addvars(summaryTableF_ODH,repmat(densityE_mean,9,1),...
    'NewVariableNames',{'ρ_e'},'After','ρ_h')
summaryNoteF_ODH = {'* mean values of ρ_k and ρ_h are given for each of the 9 combinations';...
    '* mean values of τ_s and ρ_e are given for all encounters';...
    '* predicted posterior probabilities of exploitation are evaluated for these average covariate values';...
    '* observed posterior probabiliteis are averaged for all encounters containing that combination of labels'};

summaryTableF_ODE = array2table([ceil(0.7*pairInfo.ODE{:,{'ODK','ODE'}}),...
    pairInfo.ODE{:,{'densityK_observed','densityE_observed','pE_predicted_noH','pE_predicted','pE_observed'}}],...
    'VariableNames',{'label (10^ρ_k)','label (10^ρ_e)','ρ_k','ρ_e','no history','with history','observed'});
summaryTableF_ODE(pairInfo.ODE.ODE == 200,:) = [];
summaryTableF_ODE = addvars(summaryTableF_ODE,repmat(timeOffSinceExploit_mean,9,1),...
    repmat(densityH_mean,9,1),'NewVariableNames',{'τ_s','ρ_h'},'After','ρ_k')
summaryNoteF_ODE = {'* mean values of ρ_k and ρ_e are given for each of the 9 combinations';...
    '* mean values of τ_s and ρ_h are given for all encounters';...
    '* predicted posterior probabilities of exploitation are evaluated for these average covariate values';...
    '* observed posterior probabiliteis are averaged for all encounters containing that combination of labels'};

% [~,wormList] = ismember(encounter.wormNum(indEncounter),wormNums);
% [~,conditionList] = ismember(G(indEncounter),conditionG);
% encounterList = arrayfun(@(w) 1:sum(encounter.wormNum(indEncounter)==w),unique(encounter.wormNum(indEncounter)),'UniformOutput',false);
% encounterList = [encounterList{:}]';
% conditionLabels = {'1';'1 4 7';'1 7';'1 4';'7';'4';'4 7'};
% 
% sourceTableH = table(ceil(encounter.lawnOD600(indEncounter)*0.7),wormList,encounterList,...
%     encounter.duration(indEncounter)./60,log10(encounter.duration(indEncounter)./60),...
%     encounter.velocityOn(indEncounter),log10(encounter.velocityOn(indEncounter)),...
%     encounter.exploitPosterior(indEncounter),1-encounter.exploitPosterior(indEncounter),...
%     encounter.sensePosterior(indEncounter),1-encounter.sensePosterior(indEncounter),...
%     'VariableNames',{'label','worm #','encounter #','encounter duration (min)',...
%     'log10(duration)','mean velocity on patch (μm/s)','log10(velocity)',...
%     'p(exploit)','p(explore)','p(sense)','p(non-sense)'})
% sourceNoteF = {'* posterior probability of classification as explore or exploit as estimated by the Gaussian Mixture Model described in methods and Figure 2 - supplement 6';...
%     '* posterior probability of classification as sense or non-sense as estimated by semi-supervised Quadratic Discriminant Analysis (QDA) described in methods, Figure 2I, Figure 2 - supplement 7, and Video 5';...
%     '* RGB value calculation: p(exploit)*[0 114 178] + p(explore)*p(sense)*[0 158 115] + p(explore)*p(non-sense)*[230 159 0]'};

[wormID,~] = unique(allVars(:,{'conditionNum','wormNum','encounterNum','repNum'}),'stable');
[~,wormList] = ismember(allVars.wormNum,wormNums);

sourceTableF = table(wormList,allVars.repNum,allVars.encounterNum,....
    ceil(0.7*allVars.ODK),ceil(0.7*allVars.ODH),ceil(0.7*allVars.ODE),...
    allVars.densityK,allVars.timeOffSinceExploit,allVars.densityH,allVars.densityE,...
    allVars.exploitK_GLM3,allVars.exploitK_GLM5,allVars.exploitK,...
    'VariableNames',{'worm #','replicate #','sensed encounter #',...
    'label (10^ρ_k)','label (10^ρ_h)','label (10^ρ_e)','ρ_k','τ_s','ρ_h','ρ_e',...
    'no history','with history','observed'});
sourceTableF{allVars.ODH == 200,'label (10^ρ_h)'} = 200;
sourceTableF{allVars.ODE == 200,'label (10^ρ_e)'} = 200;
sourceNoteF = {'* observed posterior probability of classification estiamted by the Gaussian Mixture Model described in methods, Figure 2H, and Figure 2 - supplement 6';...
    '* predicted posterior probabilities of exploitation are evaluated for the models with and without the history terms ρ_h and ρ_e'};

writeSourceData([saveDir,'fig4_data3.xlsx'],{'Figure 4 - Source Data 3';...
    'Figure 4F - Observed and predicted (with and without history-dependence) posterior probabilities of exploitation'},...
    {'Coefficient values for nested GLMs';'Likelihood-Ratio Test for model with and without history';...
    'Summary Statistics for ρ_k vs. ρ_h';...
    'Summary Statistics for ρ_k vs. ρ_e';'Source Data'},...
    {betaTableF;logLikeTableF;summaryTableF_ODH;summaryTableF_ODE;sourceTableF},...
    {betaNoteF;logLikeNoteF;summaryNoteF_ODH;summaryNoteF_ODE;sourceNoteF})
