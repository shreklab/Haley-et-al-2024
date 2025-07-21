%% Load Data

expName = 'foragingMini';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
saveDir = [path,'figures\FigS2\'];
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,'defineEncounter.mat'));
thresh = readtable(fullfile(path,expName,'encounterThresholds.csv'));
head = readtable([path,'foragingMini\videos\23-08-25\data\head\2023-08-25_11-15-33_2.csv']);
midpoint = readtable([path,'foragingMini\videos\23-08-25\data\midpoint\2023-08-25_11-15-33_2.csv']);

%% Figure S2A - Behavior Frames

wormNum = 214;
frameNum = 4366;
indInfo = find(cellfun(@(w) any(w == wormNum),info.wormNum));
video = getVideoInfo(info.videoFileName{indInfo});
imwrite(read(video.vidObject,frameNum),[saveDir,'figS2a.png'])

%% Figure S2B - WormLab Tracking

% Downloaded annotated frame 4366 of video of worm in Fig 3E using WormLab

%% Figure S2C - Spline Points + Patch Detection

centerpointsFile = [path,'foragingMini\videos\23-08-25\wormlab\centerpoints\2023-08-25_11-15-33_2.csv'];
wormTrack = readmatrix(centerpointsFile,'Range',[6,1]);
wormX = wormTrack(frameNum,3:2:end);
wormY = wormTrack(frameNum,4:2:end);

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
lawnOD600 = info.lawnClosestOD600{indInfo}.*info.lawnMask{indInfo};
lawnBorder = ones(size(lawnOD600));
for i = 1:length(uniqueOD600)
    lawnBorder(lawnOD600 == uniqueOD600(i)) = colorValue(i);
end
lawnImage = uint8(255.*lawnBorder);
arenaMask = uint8(info.arenaMask{indInfo});

figure('Position',[0 0 560 560])
displayImage = flip(arenaMask.*lawnImage);
displayImage(displayImage == 0) = 34; % Keynote gray
h = imshow(displayImage);
hold on
plot(wormX,wormY,'k','LineWidth',1)
scatter(wormX(1),wormY(1),25,'k.')
scaleBarSize = 2; % mm
scale = info.scale(indInfo);
pixels = info.pixels{indInfo};
scaleBar = plot(24.25*pixels/25 - [0 scaleBarSize*scale],0.75*pixels/25 + [0 0],...
    'Color','w','LineWidth',1);
scaleText = text(24.25*pixels(1)/25 - (scaleBarSize/2)*scale,pixels(2)/25,[num2str(scaleBarSize),' mm'],...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Color','w','FontName','Arial','FontSize',8);
set(gca,'YDir','normal','XLim',[1 1024],'YLim',[1 1024],...
    'XTick',[],'YTick',[],'Units','inches')
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 2])
exportgraphics(gca,[saveDir,'figS2c.pdf'],'BackgroundColor','none','ContentType','vector')

set(h,'CData',flip(read(video.vidObject,frameNum)) - 20*uint8(info.lawnMask{indInfo}));
xlim(mean(wormX) + [-100 100])
ylim(mean(wormY) + [-100 100])
exportgraphics(gca,[saveDir,'figS2c_inset.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S2D - Head-to-Midpoint Distance when Head is at Patch Edge

frameNums = [2999, 7492, 25002];
wormX = wormTrack(frameNums,3:2:end);
wormY = wormTrack(frameNums,4:2:end);

figure('Position',[0 0 560 560])
h = imshow(displayImage);
hold on
plot(wormX',wormY','k','LineWidth',1)
plot(wormX(:,[1,13])',wormY(:,[1,13])','b','LineWidth',0.75)
scatter(wormX(:,1),wormY(:,1),50,'k.')
scatter(wormX(:,13),wormY(:,13),50,'b.')
xlim(info.lawnCenters{indInfo}(1) + [-175 175])
ylim(info.lawnCenters{indInfo}(2) + [-175 175])
set(gca,'YDir','normal','XTick',[],'YTick',[],'Units','inches')
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[1 9 1.5 1.5])

subplot(555); hold on;
histogram(headMidpointDist,0:0.01:0.6,'Normalization','probability',...
    'FaceColor',colorValue.*ones(3,1))
xline(-thresh.distMidpointEnter,'b','LineWidth',1)
text(-thresh.distMidpointEnter - 0.01,0.095,['P_{50%} = ',num2str(-thresh.distMidpointEnter)],...
    'HorizontalAlignment','right','VerticalAlignment','bottom',...
    'Color','b','FontName','Arial','FontSize',8);
ylim([0 0.12])
xlim([0 0.6])
set(gca,'YDir','normal','YTick',[],'Units','inches','FontName','Arial','FontSize',8)
xlabel('head-to-midpoint distance (mm)','FontName','Arial','FontSize',8)
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[5 9 1.25 1])

subplot(5,5,20); hold on;
h = plot(head.time/60,head.distanceLawnEdge./scale,'k','LineWidth',1);
m = plot(midpoint.time/60,midpoint.distanceLawnEdge./scale,'b','LineWidth',1);
yline(0,'k','LineWidth',0.75);
yline(thresh.distMidpointEnter,'k-.','LineWidth',0.75);
yline(thresh.distMidpointMin,'k--','LineWidth',0.75);
n = xline(arrayfun(@(f) head.time(head.frameNum == f)/60,event.enter{event.wormNum == wormNum}),'g','LineWidth',0.75);
x = xline(arrayfun(@(f) head.time(head.frameNum == f)/60,event.exit{event.wormNum == wormNum}),'r','LineWidth',0.75);
set(gca,'YDir','normal','Units','inches','FontName','Arial','FontSize',8)
ax = gca; set(gca,'Position',[1 1 5 1])
ylim([-4 1])
xlabel('time (min)','FontName','Arial','FontSize',8)
ylabel('distance to patch edge (mm)','FontName','Arial','FontSize',8)
legend([h m n(1) x(1)],{'head','midpoint','enter','exit'},'Location','northeastoutside')
exportgraphics(gcf,[saveDir,'figS2d.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S2E - Head-to-Midpoint Distance when Head is at Patch Edge

frameNums = [13190, 2132, 25264];
wormX = wormTrack(frameNums,3:2:end);
wormY = wormTrack(frameNums,4:2:end);

[patchY,patchX] = find(flip(bwmorph(info.lawnMask{indInfo},'remove')));
[~,closestPatch] = min(pdist2([wormX(:,13),wormY(:,13)],[patchX,patchY]),[],2);

figure('Position',[0 0 560 560])
h = imshow(displayImage);
hold on
plot(wormX',wormY','k','LineWidth',1)
plot([wormX(:,13),patchX(closestPatch)]',[wormY(:,13),patchY(closestPatch)]','b','LineWidth',0.75)
scatter(wormX(:,1),wormY(:,1),50,'k.')
scatter(wormX(:,13),wormY(:,13),50,'b.')
xlim(info.lawnCenters{indInfo}(1) + [-175 175])
ylim(info.lawnCenters{indInfo}(2) + [-175 175])
set(gca,'YDir','normal','XTick',[],'YTick',[],'Units','inches')
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[1 9 1.5 1.5])

subplot(555); hold on;
histogram(midpointDist,-0.6:0.05:1.2,'Normalization','probability',...
    'FaceColor',colorValue.*ones(3,1))
xline(thresh.distMidpointMin,'b','LineWidth',1)
text(thresh.distMidpointMin + 0.05,0.15,['P_{1%} = ',num2str(thresh.distMidpointMin)],...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Color','b','FontName','Arial','FontSize',8);
xlim([-0.8 1.2]); ylim([0 0.19])
set(gca,'YDir','normal','YTick',[],'Units','inches','FontName','Arial','FontSize',8)
xlabel('midpoint-to-patch distance (mm)','FontName','Arial','FontSize',8)
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[5 9 1.25 1])

subplot(5,5,20); hold on;
h = plot(head.time/60,head.distanceLawnEdge./scale,'k','LineWidth',1);
m = plot(midpoint.time/60,midpoint.distanceLawnEdge./scale,'b','LineWidth',1);
yline(0,'k','LineWidth',0.75);
yline(thresh.distMidpointEnter,'k-.','LineWidth',0.75);
yline(thresh.distMidpointMin,'k--','LineWidth',0.75);
n = xline(arrayfun(@(f) head.time(head.frameNum == f)/60,event.enter{event.wormNum == wormNum}),'g','LineWidth',0.75);
x = xline(arrayfun(@(f) head.time(head.frameNum == f)/60,event.exit{event.wormNum == wormNum}),'r','LineWidth',0.75);
set(gca,'YDir','normal','Units','inches','FontName','Arial','FontSize',8)
ax = gca; set(gca,'Position',[1 1 5 1])
xlabel('time (min)','FontName','Arial','FontSize',8)
ylabel('distance to patch edge (mm)','FontName','Arial','FontSize',8)
legend([h m n(1) x(1)],{'head','midpoint','enter','exit'},'Location','northeastoutside')
xlim([3.5 5]); ylim([-4 1])
exportgraphics(gcf,[saveDir,'figS2e.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S2F - Patch-to-Midpoint Variability

figure('Position',[0 0 700 560]);
subplot(551); hold on;
varOn = vertcat(event.variabilityOn{:});
varOff = vertcat(event.variabilityOff{:});
h1 = histogram(varOn,0:0.025:1.5);
h2 = histogram(varOff(varOff >= thresh.distVarMax),0:0.025:1.5);
h3 = histogram(varOff(varOff <= thresh.distVarMax),0:0.025:1.5);
xline(thresh.distVarMax,'b','LineWidth',1)
text(thresh.distVarMax + 0.05,400,['posterior_{50%} = ',num2str(thresh.distVarMax)],...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Color','b','FontName','Arial','FontSize',8);
xlim([0 1.5]);
legend([h1 h2 h3],{'on','off (high var)','off (low var)'},'Location','northeastoutside');
set(gca,'YDir','normal','YTick',[],'Units','inches','FontName','Arial','FontSize',8)
xlabel('std(distance from patch edge)','FontName','Arial','FontSize',8)
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[1 4 1.5 1])

subplot(5,5,20); hold on;
h = plot(head.time/60,head.distanceLawnEdge./scale,'k','LineWidth',1);
m = plot(midpoint.time/60,midpoint.distanceLawnEdge./scale,'b','LineWidth',1);
yline(0,'k','LineWidth',0.75);
yline(thresh.distMidpointEnter,'k-.','LineWidth',0.75);
yline(thresh.distMidpointMin,'k--','LineWidth',0.75);
n = xline(arrayfun(@(f) head.time(head.frameNum == f)/60,event.enter{event.wormNum == wormNum}),'g','LineWidth',0.75);
x = xline(arrayfun(@(f) head.time(head.frameNum == f)/60,event.exit{event.wormNum == wormNum}),'r','LineWidth',0.75);
set(gca,'YDir','normal','Units','inches','FontName','Arial','FontSize',8)
ax = gca; set(gca,'Position',[1 1 5 1])
xlabel('time (min)','FontName','Arial','FontSize',8)
ylabel('distance to patch edge (mm)','FontName','Arial','FontSize',8)
legend([h m n(1) x(1)],{'head','midpoint','enter','exit'},'Location','northeastoutside')
xlim([6 8.5]); ylim([-4 1])
exportgraphics(gcf,[saveDir,'figS2f.pdf'],'BackgroundColor','none','ContentType','vector')