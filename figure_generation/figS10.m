%% Load Data

path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
expName = 'foragingConcentration';
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,expName,'permutePatchLocation.mat'),'onPatchShuffled');
thresh = readtable([path,'foragingMini\encounterThresholds.csv']);
saveDir = [path,'figures\FigS10\'];

permuteColor = [204 121 167]./255; % pink

%% Get condition ids

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume','growthCondition',...
    'lawnOD600','peptone'}));
conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 0.5 & ...
    ~(GID.lawnOD600 == 0.1 & GID.growthCondition == 48));
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

%% Figure 2F - p(on patch) as function of time

oneHour = (0:3599)./60;
nReps = size(onPatchShuffled,3);
smoothSize = 180; % seconds

% Calculate on-patch (midpoint w/in 0.5 mm) 
onPatchUnshuffled = false(numWorms,length(oneHour));
for i = 1:numWorms
    indData = data.wormNum == wormNums(i);
    indInfo = info.plateNum == unique(data.plateNum(indData));
    onPatchUnshuffled(i,:) = interp1(data.timeOffset(indData),...
        double(data.distanceLawnEdge(indData) >= thresh.distMidpointEnter.*info.scale(indInfo)),...
        oneHour.*60) >= 0.5;
end

% Bootstrap sample unpermuted data to obtain confidence intervals
onPatchUnshuffledBoot = nan(nReps,length(oneHour),numGroups);
for i = 1:numGroups
    onPatchUnshuffledBoot(:,:,i) = sort(bootstrp(nReps,@(X) mean(X,'omitnan'),onPatchUnshuffled(wormGroupInd == i,:)));
end
onPatchUnshuffledCI = onPatchUnshuffledBoot(round(nReps.*[0.025,0.5,0.975]),:,:);

% Get pseudorandom patch position data
onPatchShuffledExp = onPatchShuffled(ismember(unique(data.wormNum),wormNums),:,:);
onPatchShuffledMean = nan(nReps,length(oneHour),numGroups);
for i = 1:numGroups
     onPatchShuffledMean(:,:,i) = sort(reshape(mean(onPatchShuffledExp(wormGroupInd == i,:,:),1,'omitnan'),[],nReps)',1);
end
onPatchShuffledCI = onPatchShuffledMean(round(nReps.*[0.025,0.5,0.975]),:,:);

% Build contingency table
pPatchNoOverlap = nan(length(oneHour),numGroups);
pPatchNoOverlapBH = nan(length(oneHour),numGroups);
hPatchNoOverlapBH = nan(length(oneHour),numGroups);
contingency = nan(2,2,length(oneHour),numGroups);
for j = 1:numGroups
    contingency(:,:,:,j) = reshape([sum(onPatchUnshuffled(wormGroupInd == j,:) == 0,1);sum(onPatchShuffledExp(wormGroupInd == j,:,:) == 0,[1,3]);...
        sum(onPatchUnshuffled(wormGroupInd == j,:) == 1,1);sum(onPatchShuffledExp(wormGroupInd == j,:,:) == 1,[1,3])],2,2,[]); % rows = {unshuffled;shuffled} , cols = {no,yes}
    for i = 1:length(oneHour)
        % Fisher's Exact Test
        [~,pPatchNoOverlap(i,j)] = fishertest(contingency(:,:,i,j),'Tail','left');

        % Mann-Whitney U-Test
        % pPatchNoOverlap(i) = ranksum(onPatchUnshuffled(:,i),reshape(onPatchShuffled(:,i,:),[],1));
    end

    % Benjamini-Hochberg Correction
    [hPatchNoOverlapBH(:,j),pPatchNoOverlapBH(:,j)] = benjaminiHochberg(pPatchNoOverlap(:,j),0.001);
end
hPatchNoOverlapBH = double(hPatchNoOverlapBH);
hPatchNoOverlapBH(hPatchNoOverlapBH == 0) = NaN;

for j = 1:numGroups
    figure; hold on;
    patch([oneHour,flip(oneHour)],[movmean(onPatchShuffledCI(1,:,j),smoothSize,2),...
        movmean(flip(onPatchShuffledCI(3,:,j)),smoothSize,2)],...
        min(permuteColor + (1 - permuteColor).*(colorValue(10)-0.1),1),'EdgeColor','none');
    plot(oneHour,movmean(onPatchShuffledCI(2,:,j),smoothSize,2),'Color',permuteColor,'LineWidth',0.75);
    patch([oneHour,flip(oneHour)],[movmean(onPatchUnshuffledCI(1,:,j),smoothSize,2),...
        movmean(flip(onPatchUnshuffledCI(3,:,j)),smoothSize,2)],'k','EdgeColor','none',...
        'FaceAlpha',1-colorValue(10));
    plot(oneHour,movmean(onPatchUnshuffledCI(2,:,j),smoothSize,2),'k','LineWidth',0.75);
    text(mean(oneHour(~isnan(hPatchNoOverlapBH(:,j)))),0.97,'***',...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'FontName','Arial')
    plot(oneHour,hPatchNoOverlapBH(:,j).*1.05,'k','LineWidth',1)
    xlim([0 60]); ylim([0 1.2]); xticks(0:15:60); yticks([0 1])
    xlabel('time (min)','FontSize',8); ylabel('p(on patch)','FontSize',8)
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
    ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
    xlabel('time (min)','FontSize',8); ylabel('p(on patch)','FontSize',8)
    exportgraphics(gca,fullfile(saveDir,['figS10_',...
        num2str(round(10*relativeBorder(j),1,'significant'),'%.4d'),'.pdf']),...
        'BackgroundColor','none','ContentType','vector')
end