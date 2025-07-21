%% Load Data

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
data = readtable([path,expName,'\videos\22-03-18\data\midpoint\2022-03-18_13-04-43_2.csv']);
load(fullfile(path,expName,'permutePatchLocation.mat'),'onPatchShuffled');
thresh = readtable([path,'foragingMini\encounterThresholds.csv']);
saveDir = [path,'figures\FigS4\'];

%% Figure S4A - Show permutation of patch locations

wormNum = 360;
indInfo = find(cellfun(@(w) any(w == wormNum),info.wormNum));

borderAmplitude = readtable('Z:\jhaley\foragingPaper\foragingGFP\borderAmplitude.csv');
uniqueOD600 = unique(info.lawnClosestOD600{indInfo});
if uniqueOD600 > 1e-2
    indBorder = find(strcmp(borderAmplitude.peptone,info.peptone(indInfo)) & ...
        borderAmplitude.growthTimeCondition >= info.growthCondition(indInfo) & ...
        ismember(borderAmplitude.OD600,uniqueOD600) & borderAmplitude.lawnVolume == info.lawnVolume(indInfo));
    ind10 = find(strcmp(borderAmplitude.peptone,'with') & ...
        borderAmplitude.growthTimeCondition == 0 & ...
        borderAmplitude.OD600 == 10 & borderAmplitude.lawnVolume == 0.5);
    borderValue = borderAmplitude.slope(indBorder).*60.*(info.growthCondition(indInfo)+1) + ...
        borderAmplitude.intercept(indBorder);
    borderValue10 = borderAmplitude.slope(ind10).*92 + borderAmplitude.intercept(ind10);
    borderValue = 10*borderValue./borderValue10;
else
    borderValue = 1e-2;
end
colorValue = min(1 - (log(borderValue)*0.09 + 0.4),0.8);

lawnRadii = info.lawnRadii{indInfo};
pixels = info.pixels{indInfo};
lawnMask = info.lawnMask{indInfo};
arenaMask = info.arenaMask{indInfo};
noTrack = any([isnan(data.xPosition),isnan(data.yPosition)],2);

% Get x,y position of arena mask
[xArena,yArena] = find(arenaMask);
translateSpace = [xArena - pixels(1)/2,yArena - pixels(2)/2];

% Shuffle patch locations
rng(2)
lawnMaskShuffled = false(pixels);
for j = 1:length(lawnRadii)
    % Get mask of current patch
    patchMask = info.lawnMask{indInfo}.*(info.lawnClosest{indInfo} == j);

    % Center patch, then randomly rotate and translate
    patchMaskCentered = imtranslate(patchMask,size(patchMask)./2 - info.lawnCenters{indInfo}(j,:),...
        'nearest','FillValues',0);
    patchMaskRotated = imrotate(patchMaskCentered,randi(360),'nearest','crop');
    patchMaskShuffled = imtranslate(patchMaskRotated,translateSpace(randi(length(translateSpace)),:),...
        'nearest','FillValues',0);

    % Check if translation overlaps w/ other patches and arena mask
    while sum((bwdist(patchMaskShuffled) <= -thresh.distMidpointEnter.*info.scale(indInfo)) & ...
            ((bwdist(lawnMaskShuffled) <= -thresh.distMidpointEnter.*info.scale(indInfo)) | ...
            ~arenaMask),'all') > 0
        patchMaskShuffled = imtranslate(patchMaskRotated,...
            translateSpace(randi(length(translateSpace)),:),'nearest','FillValues',0);
    end
    lawnMaskShuffled = lawnMaskShuffled | patchMaskShuffled;
end

%%
figure('Position',[0 0 560 560])
lawnBorder = ones(size(lawnMask)); lawnBorder(lawnMask) = colorValue;
lawnBounds = bwdist(lawnBorder < 1) & (bwdist(lawnBorder < 1) <= -thresh.distMidpointEnter.*info.scale(indInfo));
lawnImage = uint8(255.*(lawnBorder - colorValue*lawnBounds/2));
arenaImage = uint8(arenaMask);
displayImage = flip(arenaImage.*lawnImage);
displayImage(displayImage == 0) = 34; % Keynote gray
h = imshow(displayImage);
set(gcf,'Toolbar','none')
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],...
    'XTick',[],'YTick',[],'Units','inches')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 2])
exportgraphics(gca,[saveDir,'figS4a.pdf'],'BackgroundColor','none','ContentType','vector')

figure('Position',[0 0 560 560])
lawnBorder = ones(size(lawnMaskShuffled)); lawnBorder(lawnMaskShuffled) = colorValue;
lawnBounds = bwdist(lawnBorder < 1) & (bwdist(lawnBorder < 1) <= -thresh.distMidpointEnter.*info.scale(indInfo));
lawnImage = uint8(255.*(lawnBorder - colorValue*lawnBounds/2));
arenaImage = uint8(arenaMask);
displayImage = flip(arenaImage.*lawnImage);
displayImage(displayImage == 0) = 34; % Keynote gray
h = imshow(displayImage);
set(gcf,'Toolbar','none')
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],...
    'XTick',[],'YTick',[],'Units','inches')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 2])
exportgraphics(gca,[saveDir,'figS4b.pdf'],'BackgroundColor','none','ContentType','vector')

%%
oneHour = ((1:3600)-1)./60;
indData = data.wormNum == wormNum & ~noTrack;
indPixel = sub2ind(pixels,round(data.yPosition(indData)),...
    round(data.xPosition(indData)));

thisMask = double((bwdist(lawnMask) <= ...
    -thresh.distMidpointEnter.*info.scale(indInfo)));
thisMaskShuffled = double((bwdist(lawnMaskShuffled) <= ...
    -thresh.distMidpointEnter.*info.scale(indInfo)));
onPatchWorm = interp1(data.time(indData),...
    thisMask(indPixel),oneHour.*60) >= 0.5;
onPatchShuffledWorm = interp1(data.time(indData),...
    thisMaskShuffled(indPixel),oneHour.*60) >= 0.5;

enterWorm = find(diff([0,onPatchWorm]) == 1);
exitWorm = find(diff([onPatchWorm,0]) == -1);
enterWormShuffled = find(diff([0,onPatchShuffledWorm]) == 1);
exitWormShuffled = find(diff([onPatchShuffledWorm,0]) == -1);

figure; hold on
for i = 1:length(enterWorm)
    patch([enterWorm(i).*[1 1] exitWorm(i).*[1 1]],...
        [0 1 1 0],colorValue.*[1 1 1],'EdgeColor','none')
end
xlim([0 3600]); xticks(0:600:3600); xticklabels(0:10:60); yticks([]);
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2.5 0.75])
xlabel('time (min)','FontSize',8);
exportgraphics(gca,[saveDir,'figS4e.pdf'])

figure; hold on
for i = 1:length(enterWormShuffled)
    patch([enterWormShuffled(i).*[1 1] exitWormShuffled(i).*[1 1]],...
        [0 1 1 0],colorValue.*[1 1 1],'EdgeColor','none')
end
xlim([0 3600]); xticks(0:600:3600); xticklabels(0:10:60); yticks([]);
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2.5 0.75])
xlabel('time (min)','FontSize',8);
exportgraphics(gca,[saveDir,'figS4f.pdf'])