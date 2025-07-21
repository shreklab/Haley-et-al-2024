%   LABELENCOUNTERS calculates: 
%       1. The probability of exploitation using a 2-cluster GMM of 
%           log(duration) & log(avg. velocity)
%       2. The probability of sensing using Semi-Supervised QDA of
%           min. velocity, deceleration & max. Δ velocity
%
%   Written 3/26/2024 by Jess Haley in MATLAB R2024a.
%
%   See also GETFORAGINGINFO, ANALYZEFORAGING, ANALYZEENCOUNTERS.

%% Load all encounters

path = 'Z:\jhaley\foragingPaper';
addpath(genpath(path))
expNames = {'foragingConcentration','foragingMatching','foragingMutants',...
    'foragingSensory','foragingMini'};
warning('off','stats:cvpartition:KFoldMissingGrp')

encounterAll = table();
for i = 1:length(expNames)
    load(fullfile(path,expNames{i},'encounter.mat'),'encounter');
    encounterAll = [encounterAll;encounter];
end
encounter = encounterAll; clear encounterAll

% Adjust growthCondition of lawns that were labeled 60 (foragingMini)
encounter.growthCondition(encounter.growthCondition == 60) = 48;

%% Add border amplitude value to each encounter

borderAmplitude = readtable(fullfile(path,'foragingGFP\borderAmplitude.csv'));

encounter.borderAmplitude(:) = 0;
for i = 1:size(encounter,1)

    if encounter.lawnOD600(i) > 1e-10

        % If peptone on plate use ~contains(), else just contains()
        if strcmp(encounter.peptone{i},'with')
            ind = ~contains(borderAmplitude.peptone,'out') & ...
                borderAmplitude.growthTimeCondition == encounter.growthCondition(i) & ...
                borderAmplitude.OD600 == encounter.lawnOD600(i) & ...
                borderAmplitude.lawnVolume == encounter.lawnVolume(i);
        elseif strcmp(encounter.peptone{i},'without')
            ind = contains(borderAmplitude.peptone,'out') & ...
                borderAmplitude.growthTimeCondition == encounter.growthCondition(i) & ...
                borderAmplitude.OD600 == encounter.lawnOD600(i) & ...
                borderAmplitude.lawnVolume == encounter.lawnVolume(i);
        end

        % If no value found, use NaN
        if sum(ind) == 0
            encounter.borderAmplitude(i) = NaN;
            continue
        end

        % Get border amplitude value via interpolation of timeExp
        encounter.borderAmplitude(i) = borderAmplitude.intercept(ind) + ...
            borderAmplitude.slope(ind)*encounter.lawnGrowth(i);
    end
end

encounter.borderAmplitude(isnan(encounter.enter)) = NaN;

% Get border amplitude of growth plates
growthInd = borderAmplitude.lawnVolume == 200;
encounter.borderAmplitudeGrowth = borderAmplitude.intercept(growthInd) + ...
    borderAmplitude.slope(growthInd)*minutes(encounter.growthLawnGrowth);

%% Cluster exploit using GMM of log(duration) & log(avg. velocity)

% Get cluster variables and group by experiment name and concentration
clusterVars = [log10(encounter.duration./60),log10(encounter.velocityOn)];
clusterGroups = findgroups(encounter(:,{'expName','lawnVolume','lawnOD600','growthCondition'}));
include = sum(isnan(clusterVars),2) == 0 & ~encounter.exclude;
clusterVarsInclude = clusterVars(include,:);
clusterGroups = clusterGroups(include,:);

% Define cluster colors
exploreColor = [0 158 115]./255; % green
exploitColor = [0 114 178]./255; % blue

% Define parameters
nCluster = 2;
nReps = 1000; % # repetitions
k = 5; % k-fold cross-validation
c = cvpartition(clusterGroups,'KFold',k);
alpha = 0:0.005:0.05; % regularization values

% Initialize variables
muExploit = nan(nCluster,nCluster,nReps,length(alpha),k);
sigmaExploit = nan(nCluster,nCluster,nCluster,nReps,length(alpha),k);
proportionExploit = nan(nCluster,nReps,length(alpha),k);
posteriorVarExploit = nan(nReps,length(alpha),k);

for n = 1:nReps
    % Repartition
    c = repartition(c);

    for i = 1:length(alpha)
        for j = 1:k

            % Fit Gaussian Mixture Model to k-fold training set
            clusterModel = fitgmdist(clusterVarsInclude(training(c,j),:),nCluster,...
                'RegularizationValue',alpha(i),...
                'Replicates',1,'Options',statset('MaxIter',10000));

            % Save parameters - make sure cluster 2 is exploit
            if diff(clusterModel.mu(:,1)) > 0
                muExploit(:,:,n,i,j) = clusterModel.mu;
                sigmaExploit(:,:,:,n,i,j) = clusterModel.Sigma;
                proportionExploit(:,n,i,j) = clusterModel.ComponentProportion;
            else
                muExploit(:,:,n,i,j) = flip(clusterModel.mu,1);
                sigmaExploit(:,:,:,n,i,j) = flip(clusterModel.Sigma,3);
                proportionExploit(:,n,i,j) = flip(clusterModel.ComponentProportion);
            end

            % Calculate posterior variance using k-fold test set
            posteriorProbExploit = posterior(clusterModel,clusterVarsInclude(test(c,j),:));
            posteriorVarExploit(n,i,j) = mean(prod(posteriorProbExploit,2));
        end
    end
end

% Plot posterior variance as a function of alpha
figure; hold on;
for n = 1:nReps
    scatter(alpha,permute(posteriorVarExploit(n,:,:),[3 2 1]),'filled')
end
scatter(alpha,mean(posteriorVarExploit,[1,3]),'k','filled')
xlabel('α'); ylabel('posterior variance')

% Get parameters of the best model which has the lowest variance
[~,bestModel] = min(mean(posteriorVarExploit,[1,3]));
bestMu = mean(muExploit(:,:,:,bestModel,:),3:5);
bestSigma = mean(sigmaExploit(:,:,:,:,bestModel,:),4:6);
bestProportion = mean(proportionExploit(:,:,bestModel,:),2:4);

% Create best model
exploitModel = gmdistribution(bestMu,bestSigma,bestProportion);
posteriorProbExploit = posterior(exploitModel,clusterVars);
[~,exploreCluster] = min(bestMu(:,1));
[~,exploitCluster] = max(bestMu(:,1));

% Get individual gaussians for explore and exploit
exploitDist = gmdistribution(bestMu(exploitCluster,:),bestSigma(:,:,exploitCluster));
exploreDist = gmdistribution(bestMu(exploreCluster,:),bestSigma(:,:,exploreCluster));

% Plot posterior probability for best model
figure('Position',[400 400 520 520]); hold on;
[X,Y] = meshgrid(linspace(-2,3,1000),linspace(1.2,2.8,1000));
Z = reshape(pdf(exploitDist,[X(:),Y(:)]),size(X));
contourf(X,Y,Z,0.7*max(Z,[],'all').*10.^(-2:0),'FaceAlpha',0.2,'EdgeAlpha',0);
Z = reshape(pdf(exploreDist,[X(:),Y(:)]),size(X));
contourf(X,Y,Z,0.7*max(Z,[],'all').*10.^(-2:0),'FaceAlpha',0.2,'EdgeAlpha',0);
colormap(flipud(gray)); set(gca,'ColorScale','log');
scatter(clusterVars(:,1),clusterVars(:,2),5,'CData',...
    posteriorProbExploit(:,exploreCluster).*exploreColor + posteriorProbExploit(:,exploitCluster).*exploitColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none')
title(['α = ',num2str(alpha(bestModel))])
xlabel('log(time on patch)')
ylabel('log(velocity on patch)')

% Assign posterior probability of exploiting to each encounter
encounter.exploitPosterior = posteriorProbExploit(:,exploitCluster);

%% Cluster sensing using Semi-Supervised QDA of min. velocity, deceleration & max. Δ velocity

% Get cluster variables and group by experiment name and concentration
clusterVars = [encounter.velocityOnMin,...
    encounter.velocityBeforeEnter-encounter.velocityOnMin,encounter.decelerate];
[clusterGroups,GID] = findgroups(encounter(:,...
    {'expName','peptone','lawnVolume','lawnOD600','growthCondition'}));
nonsenseGroup = find(GID.lawnOD600 == 1e-10);
startsOn = encounter.timeEnter <= 10; % exclude encounters where velocityBeforeEnter is censored
minVelocityOnly = ~isnan(clusterVars(:,1)) & ...
    (any(isnan(clusterVars(:,2:3)),2) | startsOn); % exclude encounters where only clusterVar1 is reliable
include = ~any(isnan(clusterVars),2) & ~encounter.exclude & ~minVelocityOnly;

% Define cluster colors
nonSenseColor = [230 159 0]./255; % orange
exploreColor = [0 158 115]./255; % green

% Define parameters
nCluster = 2;
nReps = 1000; % # repetitions
subSample = 1.*(1:-0.1:0.1); % fraction of data to subsample

% Subsample indices for inclusion
getSubset = @(indices,sampleFrac) indices & (rand(size(indices)) < sampleFrac);

% Initialize variables
muSense = nan(nCluster,size(clusterVars,2),nReps,length(subSample));
sigmaSense = nan(size(clusterVars,2),size(clusterVars,2),nCluster,nReps,length(subSample));
priorSense = nan(nCluster,nReps,length(subSample));
coeffConstant = nan(nReps,length(subSample));
coeffLinear =  nan(size(clusterVars,2),nReps,length(subSample));
coeffQuadratic = nan(size(clusterVars,2),size(clusterVars,2),nReps,length(subSample));
posteriorVarSense = nan(nReps,length(subSample));
posteriorChiSquare = nan(nReps,length(subSample));

for n = 1:nReps
    
    % Use posterior probability of exploit to construct sense training set
    indSense = rand(height(encounter),1) <= encounter.exploitPosterior;
    indNonsense = ismember(clusterGroups,nonsenseGroup);
    trainSet = include & (indSense | indNonsense);

        for j = 1:length(subSample)

            % Subsample training set to test sensitivity of model
            trainSub = getSubset(trainSet,subSample(j));

            % Fit semi-self supervised QDA to training set
            discrimModel = fitsemiself(clusterVars(trainSub,:),...
                indNonsense(trainSub),clusterVars(~trainSub & include,:),'Learner',...
                templateDiscriminant('DiscrimType','quadratic','Cost',[0 1;1 0]));

            % Save parameters - make sure cluster 1 is sense (i.e. lower min velocity)
            if diff(discrimModel.Learner.Mu(:,1)) > 0
                muSense(:,:,n,j) = discrimModel.Learner.Mu;
                sigmaSense(:,:,:,n,j) = discrimModel.Learner.Sigma;
                priorSense(:,n,j) = discrimModel.Learner.Prior;
                coeffConstant(n,j) = discrimModel.Learner.Coeffs(2,1).Const;
                coeffLinear(:,n,j) = discrimModel.Learner.Coeffs(2,1).Linear; 
                coeffQuadratic(:,:,n,j) = discrimModel.Learner.Coeffs(2,1).Quadratic;
                senseCluster = 2;
            else
                muSense(:,:,n,j) = flip(discrimModel.Learner.Mu,1);
                sigmaSense(:,:,:,n,j) = flip(discrimModel.Learner.Sigma,3);
                priorSense(:,n,j) = flip(discrimModel.Learner.Prior,1);
                coeffConstant(n,j) = discrimModel.Learner.Coeffs(1,2).Const;
                coeffLinear(:,n,j) = discrimModel.Learner.Coeffs(1,2).Linear; 
                coeffQuadratic(:,:,n,j) = discrimModel.Learner.Coeffs(1,2).Quadratic;
                senseCluster = 1;
            end

            % Calculate posterior variance of test set
            [~,posteriorProbSense] = predict(discrimModel,clusterVars(~trainSub & include,:));
            posteriorVarSense(n,j) = mean(prod(posteriorProbSense,2));

            % Calculate variance of test set posteriors b/w subsamples and
            % full-sized training set
            [~,posteriorProbSense] = predict(discrimModel,clusterVars(~trainSet & include,:));
            if subSample(j) == 1
                posteriorProbFull = posteriorProbSense(:,senseCluster);
            end
            posteriorChiSquare(n,j) = sum((posteriorProbSense(:,senseCluster) - posteriorProbFull).^2./...
                posteriorProbFull);
        end
    
end

% Plot posterior variance as a function of subsample size
figure;
subplot(211); plot(subSample,posteriorVarSense)
set(gca,'XDir','reverse'); xlabel('fraction of training set'); ylabel('posterior variance');
subplot(212); plot(subSample,posteriorChiSquare)
set(gca,'XDir','reverse'); xlabel('fraction of training set'); ylabel('posterior chi-squared');

% Plot variance of quadratic coefficients as a function of subsample size
figure;
subplot(521); plot(subSample,coeffConstant); title('β_0'); set(gca,'XDir','reverse');
for i = 1:size(coeffLinear,1)
    subplot(5,2,i+1);
    plot(subSample,squeeze(coeffLinear(i,:,:)))
    title(['β*x_',num2str(i)]); set(gca,'XDir','reverse');
end
for i = 1:size(coeffQuadratic,1)
    subplot(5,2,i+4);
    plot(subSample,squeeze(coeffQuadratic(i,i,:,:)))
    title(['β*x_',num2str(i),'^2']); set(gca,'XDir','reverse');
end
subplot(5,2,8); plot(subSample,...
    squeeze(coeffQuadratic(1,2,:,:) + coeffQuadratic(2,1,:,:)))
title('β*x_1*x_2'); set(gca,'XDir','reverse');
subplot(5,2,9); plot(subSample,...
    squeeze(coeffQuadratic(1,3,:,:) + coeffQuadratic(3,1,:,:)))
title('β*x_1*x_3'); set(gca,'XDir','reverse');
subplot(5,2,10); plot(subSample,...
    squeeze(coeffQuadratic(2,3,:,:) + coeffQuadratic(3,2,:,:)))
title('β*x_2*x_3'); set(gca,'XDir','reverse');

% Get parameters of the full-size training set model
fullModel = 1;
bestMu = mean(muSense(:,:,:,fullModel),3:4);
bestSigma = mean(sigmaSense(:,:,:,:,fullModel),4:5);
bestPrior = mean(priorSense(:,:,fullModel),2:3);

% Create best model
senseModel = makecdiscr(bestMu,bestSigma,'Prior',bestPrior);
[~,posteriorProbSense] = predict(senseModel,clusterVars);
[~,senseCluster] = min(bestMu(:,1));
[~,nonsenseCluster] = max(bestMu(:,1));

% Plot posterior probability for best model
figure; subplot(211);
scatter3(clusterVars(trainSet,1),clusterVars(trainSet,2),clusterVars(trainSet,3),...
    5,'CData',indSense(trainSet).*[1 1 1],...
    'MarkerFaceColor','flat','MarkerEdgeColor','k','MarkerFaceAlpha',1)
view(0,90); title('labeled training set')
subplot(212); hold on
scatter3(clusterVars(:,1),clusterVars(:,2),clusterVars(:,3),5,'CData',...
    posteriorProbSense(:,senseCluster).*exploreColor + posteriorProbSense(:,nonsenseCluster).*nonSenseColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
K = senseModel.Coeffs(1,2).Const;
L = senseModel.Coeffs(1,2).Linear; 
Q = senseModel.Coeffs(1,2).Quadratic;
f = @(x1,x2,x3) K + L(1)*x1 + L(2)*x2 + L(3)*x3 + ...
    Q(1,1)*x1.^2 + Q(2,2)*x2.^2 + Q(3,3)*x3.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + (Q(1,3)+Q(3,1))*x1.*x3 + (Q(2,3)+Q(3,2))*x2.*x3;
fimplicit3(f,[0 500 0 500 -60 40],'EdgeColor','none','FaceAlpha',0.5,'FaceColor','k'); legend('off')
view(0,90); title('posterior probabilities')

%% Estimate posterior probabilities from min. velocity only via marginalization

% Get min velocities and estimate posterior
clusterVar1 = clusterVars(minVelocityOnly,1);
posteriorProbSenseMarginal = nan(size(clusterVar1));
for i = 1:length(clusterVar1)
    posteriorProbSenseMarginal(i) = estimatePosterior(senseModel,clusterVars(include,:),clusterVar1(i));
end

% Plot change in posterior probabilities
figure; hold on
gscatter(posteriorProbSense(minVelocityOnly,senseCluster),...
    posteriorProbSenseMarginal,...
    encounter.timeEnter(minVelocityOnly) > 0)
plot([0 1],[0 1],'k--')
xlabel('censored p(sense)'); ylabel('marginalized p(sense)')
legend('time enter == 0','time enter <= 10'); xlim([0 1])

% Plot posterior probability for best model
figure('Position',[400 400 520 520]); hold on;
scatter(clusterVars(:,1),clusterVars(:,2),15,'CData',...
    posteriorProbSense(:,senseCluster).*exploreColor + posteriorProbSense(:,nonsenseCluster).*nonSenseColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',0.1)
% scatter(clusterVars(minVelocityOnly,1),clusterVars(minVelocityOnly,2),10,'CData',...
%     posteriorProbSenseMarginal.*exploreColor + (1-posteriorProbSenseMarginal).*nonSenseColor,...
%     'MarkerFaceColor','flat','MarkerEdgeColor','k')
xlabel('min. velocity on patch (μm/s)')
ylabel('max. Δ velocity (μm/s)')

% Assign posterior probability of sensing to each encounter
encounter.sensePosterior = posteriorProbSense(:,senseCluster);
encounter.sensePosterior(minVelocityOnly) = posteriorProbSenseMarginal;

%% Compute probabilities

% Marginal Probabilities
marginal.exploit = mean(encounter.exploitPosterior(~encounter.exclude),'omitnan');
marginal.explore = 1 - marginal.exploit;
marginal.sense = mean(encounter.sensePosterior(~encounter.exclude),'omitnan');
marginal.nonsense = 1 - marginal.sense;

% Joint Probabilities
joint.sense_exploit = mean(encounter.exploitPosterior(~encounter.exclude).*...
    encounter.sensePosterior(~encounter.exclude),'omitnan');
joint.nonsense_exploit = mean(encounter.exploitPosterior(~encounter.exclude).*...
    (1 - encounter.sensePosterior(~encounter.exclude)),'omitnan');
joint.sense_explore = mean((1 - encounter.exploitPosterior(~encounter.exclude)).*...
    encounter.sensePosterior(~encounter.exclude),'omitnan');
joint.nonsense_explore = mean((1 - encounter.exploitPosterior(~encounter.exclude)).*...
    (1 - encounter.sensePosterior(~encounter.exclude)),'omitnan');

% Conditional Probabilities
conditional.sense_exploit = joint.sense_exploit/marginal.exploit;
conditional.nonsense_exploit = joint.nonsense_exploit/marginal.exploit;
conditional.sense_explore = joint.sense_explore/marginal.explore;
conditional.nonsense_explore = joint.nonsense_explore/marginal.explore;
conditional.exploit_sense = joint.sense_exploit/marginal.sense;
conditional.exploit_nonsense = joint.nonsense_exploit/marginal.nonsense;
conditional.explore_sense = joint.sense_explore/marginal.sense;
conditional.explore_nonsense = joint.nonsense_explore/marginal.nonsense;

%% Label encounters via hard clustering

encounter.label(:) = {''};
encounter.label(encounter.exploitPosterior >= 0.5) = {'exploit'};
encounter.label(encounter.sensePosterior >= 0.5 & encounter.exploitPosterior < 0.5) = {'sample'};
encounter.label(encounter.sensePosterior < 0.5 & encounter.exploitPosterior < 0.5) = {'searchOn'};
encounter.label(encounter.sensePosterior < 0.05 & encounter.exploitPosterior < 0.5 & ...
    encounter.distanceOnMax < 0) = {'searchOff'};

%% Save analyses

save(fullfile(path,'encounter.mat'),'encounter','-v7.3');
save(fullfile(path,'labelEncounters.mat'));