%% Load all info
path = 'Z:\jhaley\foragingPaper';
addpath(genpath(path))
expNames = {'foragingConcentration','foragingMatching','foragingMutants',...
    'foragingSensory','foragingMini'};
load(fullfile(path,'encounter.mat'),'encounter');
%%
infoAll = table();
for i = 1:length(expNames)
    load(fullfile(path,expNames{i},'experimentInfo.mat'),'info');
    info.expName(:) = expNames(i);
    if width(info) ~= width(infoAll) && width(infoAll) > 0
        % Add to infoAll
        missingVars = info.Properties.VariableNames(~ismember(info.Properties.VariableNames,infoAll.Properties.VariableNames));
        for v = 1:length(missingVars)
            if iscell(info.(missingVars{v}))
                infoAll.(missingVars{v})(:) = {''};
            elseif isduration(info.(missingVars{v}))
                infoAll.(missingVars{v})(:) = duration(0,0,0);
            elseif isdatetime(info.(missingVars{v}))
                infoAll.(missingVars{v})(:) = NaT;
            else
                infoAll.(missingVars{v})(:) = NaN;
            end
        end

        % Add to info
        missingVars = infoAll.Properties.VariableNames(~ismember(infoAll.Properties.VariableNames,info.Properties.VariableNames));
        for v = 1:length(missingVars)
            if iscell(infoAll.(missingVars{v}))
                info.(missingVars{v})(:) = {''};
            elseif isduration(infoAll.(missingVars{v}))
                info.(missingVars{v})(:) = duration(0,0,0);
            elseif isdatetime(infoAll.(missingVars{v}))
                info.(missingVars{v})(:) = NaT;
            else
                info.(missingVars{v})(:) = NaN;
            end
        end
    end
    if ~iscell(info.OD600)
        info.OD600 = num2cell(info.OD600);
    end
    infoAll = [infoAll;info];
end
info = infoAll; clear infoAll missingVars

%% Patch size

diameterSmall = mean(vertcat(info.lawnDiameterMean{~info.exclude & info.lawnVolume == 0.5}))
diameterLarge = mean(vertcat(info.lawnDiameterMean{~info.exclude & info.lawnVolume == 20}))

%% Bacterial Cultures - mean OD600, CFU, and plate storage

OD600 = [mean(info.OD600Real(~info.exclude)), std(info.OD600Real(~info.exclude))] % measured OD = 10
CFU = 2*10^7*[mean(info.CFU(~info.exclude),'omitnan'),...
    std(info.CFU(~info.exclude),'omitnan')] % colonies in 1 mL of OD = 10 solution
plateStorage = days(info.timeRecord(~info.exclude) - info.timeColdRoom(~info.exclude));
plateStorage = [mean(plateStorage,'omitnan'), std(plateStorage,'omitnan')] % days in cold room

%% Nematode Cultures - experiment dates, nematode growth

earliestExp = min(info.growthTimePicked(~info.exclude))
latestExp = max(info.timeRecord(~info.exclude))
growthTime = [mean(hours(info.wormGrowth(~info.exclude)),'omitnan'),...
    std(hours(info.wormGrowth(~info.exclude)),'omitnan')]

%% Assay Preparation - bacterial growth

ind = ~info.exclude & strcmp(info.expName,'foragingConcentration') & ...
    strcmp(info.condition,'grid') & info.lawnVolume == 0.5 & ...
        ~(strcmp(info.OD600Label,'0.10') & info.growthCondition == 48);
infoConcentration = info(ind,:);
infoConcentration.OD600 = cellfun(@str2num,infoConcentration.OD600Label);
[infoConcentration.conditionNum,GID] = findgroups(infoConcentration(:,{'growthCondition','OD600'}));

metrics = {'lawnGrowth','wormGrowth','OD600Real','CFU','temp','humidity'};
aov = cell(size(metrics));
for m = 1:length(metrics)
    GID.(metrics{m}) = splitapply(@nanmean,info.(metrics{m})(ind),G);
    if isduration(infoConcentration.(metrics{m}))
        infoConcentration.(metrics{m}) = hours(infoConcentration.(metrics{m}));
    end
    aov{m} = anova(infoConcentration,metrics{m},'FactorNames','conditionNum');
end
GID

temp = [mean(info.temp(~info.exclude),'omitnan'),std(info.temp(~info.exclude),'omitnan')] % celsius
humidity = [mean(info.humidity(~info.exclude),'omitnan'),std(info.humidity(~info.exclude),'omitnan')] % percent
bacterialGrowth_00 = hours(info.lawnGrowth(~info.exclude & info.growthCondition == 0));
bacterialGrowth_00 = [mean(bacterialGrowth_00,'omitnan'),std(bacterialGrowth_00,'omitnan')] % hours
bacterialGrowth_12 = hours(info.lawnGrowth(~info.exclude & info.growthCondition == 12));
bacterialGrowth_12 = [mean(bacterialGrowth_12,'omitnan'),std(bacterialGrowth_12,'omitnan')] % hours
bacterialGrowth_48 = hours(info.lawnGrowth(~info.exclude & info.growthCondition == 48));
bacterialGrowth_48 = [mean(bacterialGrowth_48,'omitnan'),std(bacterialGrowth_48,'omitnan')] % hours

%% Resolution

spatialRes = [mean(info.scale(~info.exclude & ~strcmp(info.expName,'foragingMini'))),...
    std(info.scale(~info.exclude & ~strcmp(info.expName,'foragingMini')))]
spatialRes_mini = [mean(info.scale(~info.exclude & strcmp(info.expName,'foragingMini'))),...
    std(info.scale(~info.exclude & 3772strcmp(info.expName,'foragingMini')))]

%% Sense-nonsense labeling

load(fullfile(path,'labelEncounters.mat'),'clusterGroups','GID','nonsenseGroup',...
    'minVelocityOnly','senseCluster','posteriorProbSense','posteriorProbSenseMarginal');

nReps = 1000;
indSense = rand(height(encounter),nReps) <= encounter.exploitPosterior & ~encounter.exclude;
nonsensePosterior = encounter.sensePosterior(ismember(clusterGroups,nonsenseGroup));
sensePosterior = repmat(encounter.sensePosterior,1,nReps);
sensePosterior = sensePosterior(indSense);
senseFalsePos = mean(nonsensePosterior)
senseFalseNeg = 1 - mean(sensePosterior)

clusterVars = [encounter.velocityOnMin,...
    encounter.velocityBeforeEnter-encounter.velocityOnMin,encounter.decelerate];
numCensored = [sum(minVelocityOnly), sum(~isnan(clusterVars(:,1)))]

posteriorDelta = posteriorProbSenseMarginal - posteriorProbSense(minVelocityOnly,senseCluster);
[mean(posteriorDelta,'omitnan'), std(posteriorDelta,'omitnan')]