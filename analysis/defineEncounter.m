%   DEFINEENCOUNTER uses high-resolution body segment locations of animals
%   exploring environments with a single small patch of bacteria to define
%   an encounter event for use when only the mid-point of the animal can be
%   reliably tracked.
%
%   Written 3/8/2024 by Jess Haley in MATLAB R2023b.
%
%   See also GETFORAGINGINFO, ANALYZEWORMLABTRACKS, ANALYZEENCOUNTERS.

%% Load Data

path = 'Z:\jhaley\foragingPaper\foragingMini\';
load([path,'experimentInfo.mat'],'info');
load([path,'head.mat'],'data'); head = data;
load([path,'midpoint.mat'],'data'); midpoint = data;
load([path,'tail.mat'],'data'); tail = data;
clear data
wormNums = unique(head.wormNum);
warning('off','MATLAB:table:RowsAddedExistingVars')

%% Midpoint Distance Thresholds

dataScale = arrayfun(@(p) repmat(info.scale(info.plateNum == p),...
    sum(head.plateNum == p),1),unique(head.plateNum),'UniformOutput',false);
scale = vertcat(dataScale{:});

% Calculate average distance from head to midpoint when the head is
% touching the patch edge
ind = find(abs(head.distanceLawnEdge) <= 1);
headMidpointDist = sqrt((head.xPosition(ind) - midpoint.xPosition(ind)).^2 + ...
    (head.yPosition(ind) - midpoint.yPosition(ind)).^2)./scale(ind);
midpointThresh1 = prctile(headMidpointDist,50); % mm

% Calculate minimum distance midpoint is from patch edge when head is on
ind = find(head.distanceLawnEdge >= 0 & head.closestOD600 >= 1);
midpointDist = midpoint.distanceLawnEdge(ind)./scale(ind);
midpointThresh2 = prctile(midpointDist,1); % mm

figure('Position',[0 0 560 420]);
subplot(211); hold on
histogram(headMidpointDist,'Normalization','probability')
xline(midpointThresh1,'r'); xlim([0 0.7])
xlabel('Distance from Head-to-Midpoint')

subplot(212); hold on
histogram(midpointDist,'Normalization','probability')
xline(midpointThresh2,'r');
xlabel('Distance from Midpoint to Patch Edge')

%% Distance Variability

event = table(); 
for i = 1:length(wormNums)
    % Get indices and worm info
    ind = midpoint.wormNum == wormNums(i);
    indInfo = find(cellfun(@(w) ismember(w,wormNums(i)),info.wormNum));
    event.wormNum(i) = wormNums(i);
    event.OD600(i) = info.OD600(indInfo);
    event.growthCondition(i) = info.growthCondition(indInfo);
    
    % Find enter and exit events based on midpoint
    % enterEvent = find(diff([0;head.distanceLawnEdge(ind)./info.scale(indInfo) >= 0])==1);
    % exitEvent = find(diff([head.distanceLawnEdge(ind)./info.scale(indInfo) >= 0;0])==-1);
    enterEvent = find(diff([0;midpoint.distanceLawnEdge(ind)./info.scale(indInfo) >= -midpointThresh1])==1);
    exitEvent = find(diff([midpoint.distanceLawnEdge(ind)./info.scale(indInfo) >= -midpointThresh1;0])==-1);

    % Calculate the standard deviation of the distance from patch edge
    % during on and off patch events
    distLawnEdge = midpoint.distanceLawnEdge(ind);
    timeOffset = midpoint.timeOffset(ind);
    event.enter{i} = enterEvent;
    event.exit{i} = exitEvent;
    event.variabilityOn{i} = arrayfun(@(enter,exit) std(distLawnEdge(enter:exit)./...
        info.scale(indInfo)),enterEvent,exitEvent);
    event.variabilityOff{i} = arrayfun(@(enter,exit) std(distLawnEdge(exit:enter)./...
        info.scale(indInfo)),enterEvent(2:end),exitEvent(1:end-1));
    event.midpointOff{i} = arrayfun(@(enter,exit) max(distLawnEdge(exit:enter)./...
        info.scale(indInfo)),enterEvent(2:end),exitEvent(1:end-1));
    event.durationBeforeOff{i} = timeOffset(exitEvent(1:end-1))-timeOffset(enterEvent(1:end-1));
    event.durationOff{i} = timeOffset(enterEvent(2:end))-timeOffset(exitEvent(1:end-1));
    event.wormNumOff{i} = wormNums(i).*ones(length(exitEvent)-1,1);
end

figure('Position',[0 0 560 420]); hold on
histogram(vertcat(event.variabilityOn{:}),0:0.025:1.5)
histogram(vertcat(event.variabilityOff{:}),0:0.025:1.5)
legend({'on patch','off patch'}); xlabel('std(distance from patch edge)');

figure('Position',[0 0 560 420]);
scatter(vertcat(event.variabilityOff{:}),vertcat(event.durationBeforeOff{:}))
yscale('log'); xlabel('std(distance from patch edge)'); ylabel('log(time on previous patch)')

figure('Position',[0 0 560 420]);
gscatter(vertcat(event.variabilityOff{:}),vertcat(event.durationOff{:}),vertcat(event.wormNumOff{:}))
xlabel('std(distance from patch edge)'); ylabel('log(time off patch)'); legend('off')

% Get cluster variables and group by experiment name and concentration
clusterVars = vertcat(event.variabilityOff{:});

% Define parameters
nCluster = 2;
nReps = 1000; % # repetitions
k = 5; % k-fold cross-validation
c = cvpartition(length(clusterVars),'KFold',k);
alpha = 0;%:0.001:0.05; % regularization values

% Initialize variables
mu = nan(nCluster,size(clusterVars,2),nReps,length(alpha),k);
sigma = nan(size(clusterVars,2),size(clusterVars,2),nCluster,nReps,length(alpha),k);
proportion = nan(nCluster,nReps,length(alpha),k);
posteriorVar = nan(nReps,length(alpha),k);

for n = 1:nReps
    % Repartition
    c = repartition(c);

    for i = 1:length(alpha)
        for j = 1:k

            % Fit Gaussian Mixture Model to k-fold training set
            clusterModel = fitgmdist(clusterVars(training(c,j),:),nCluster,...
                'RegularizationValue',alpha(i),...
                'Replicates',1,'Options',statset('MaxIter',10000));

            % Save parameters - make sure cluster 2 is larger
            if diff(clusterModel.mu(:,1)) > 0
                mu(:,:,n,i,j) = clusterModel.mu;
                sigma(:,:,:,n,i,j) = clusterModel.Sigma;
                proportion(:,n,i,j) = clusterModel.ComponentProportion;
            else
                mu(:,:,n,i,j) = flip(clusterModel.mu,1);
                sigma(:,:,:,n,i,j) = flip(clusterModel.Sigma,3);
                proportion(:,n,i,j) = flip(clusterModel.ComponentProportion);
            end

            % Calculate posterior variance using k-fold test set
            posteriorProb = posterior(clusterModel,clusterVars(test(c,j),:));
            posteriorVar(n,i,j) = mean(prod(posteriorProb,2));
        end
    end
end

% Plot posterior variance as a function of alpha
figure; hold on;
for n = 1:nReps
    scatter(alpha,permute(posteriorVar(n,:,:),[3 2 1]),'filled')
end
scatter(alpha,mean(posteriorVar,[1,3]),'k','filled')
xlabel('α'); ylabel('posterior variance')

% Get parameters of the best model which has the lowest variance
[~,bestModel] = min(mean(posteriorVar,[1,3]));
bestAlpha = alpha(bestModel);
bestMu = mean(mu(:,:,:,bestModel,:),3:5);
bestSigma = mean(sigma(:,:,:,:,bestModel,:),4:6);
bestProportion = mean(proportion(:,:,bestModel,:),2:4);

% Create best model
GMModel = gmdistribution(bestMu,bestSigma,bestProportion);
posteriorProb = posterior(GMModel,clusterVars);
[~,lowVarCluster] = min(bestMu(:,1));
[~,indCluster] = min(posteriorProb,[],2);

figure('Position',[0 0 560 420]); hold on
histogram(clusterVars(posteriorProb(:,lowVarCluster) >= 0.5,1),0:0.025:1.5)
histogram(clusterVars(posteriorProb(:,lowVarCluster) < 0.5,1),0:0.025:1.5)
xlabel('var(distance from patch edge)');

[posteriorProb,ind] = unique(posteriorProb(:,lowVarCluster));
variabilityThresh = interp1(posteriorProb,clusterVars(ind,1),0.5);
xline(variabilityThresh,'k'); legend({'low variability','high variability'}); 

%% Export Results

thresholds = table(-midpointThresh1,midpointThresh2,variabilityThresh,...
    'VariableNames',{'distMidpointEnter','distMidpointMin','distVarMax'})
writetable(thresholds,[path,'encounterThresholds.csv']);
save([path,'defineEncounter.mat'],'headMidpointDist','midpointDist','event','GMModel','posteriorProb','-v7.3');
