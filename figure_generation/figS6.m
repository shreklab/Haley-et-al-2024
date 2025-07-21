% Load data

path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'foragingGFP','analyzeGFP_24-04-04.mat'),...
    'background','backgroundInfo','data','info',...
    'lawnAnalysis','lawnProfiles','metaData','template');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
saveDir = [path,'figures\FigS6\'];

% Exclude image(s) and plate(s) from analysis
excludeImages = {'20231230_101739.czi'};
excludePlates = [64];
excludeAnalysis = [cellfun(@(f) strcmp(lawnAnalysis.fileName,f),...
    excludeImages,'UniformOutput',false),...
    arrayfun(@(p) lawnAnalysis.plateNum == p,excludePlates,...
    'UniformOutput',false)];
lawnAnalysis.include = sum([excludeAnalysis{:}],2) == 0;

% Exclude outliers due to dust, poor lawn detection, etc.
lawnAnalysis.include((lawnAnalysis.circularity < 0.75 & ...
    lawnAnalysis.lawnVolume == 0.5) | ...
    lawnAnalysis.borderCenterRatio < 1 | ...
    lawnAnalysis.borderCenterRatio > 35) = 0;

% Exclude patches with poor outer edge detection
imageNums = unique(lawnAnalysis.imageNum);
for i = 1:length(imageNums)
    while diff(prctile(lawnAnalysis.xPeak(lawnAnalysis.imageNum == imageNums(i) & ...
            lawnAnalysis.include),[0 100])) > 0.1
        ind = find(lawnAnalysis.imageNum == imageNums(i) & lawnAnalysis.include);
        [~,worstFit] = min(lawnAnalysis.xPeak(ind));
        lawnAnalysis.include(ind(worstFit)) = 0;
    end
end

% Get conditions
[lawnAnalysis.condition,conditions] = ...
    findgroups(lawnAnalysis(:,{'peptone','lawnVolume','growthTimeCondition',...
    'OD600'}));
condToAnalyze = find(ismember(conditions.growthTimeCondition,[0,48]) | ...
    (conditions.growthTimeCondition == 12 & conditions.OD600 == 1 & conditions.lawnVolume == 0.5));

%% Figure S6A - Example patch image

% OD600 = 10, expNum = 3, expTime = 2 hr
imageNum = 749;
imageNumBack = metaData.backgroundImageNum(imageNum);
imageNumBright = metaData.brightfieldImageNum(imageNum);
scale = metaData.scale(imageNum);

% Raw brightfield image
figure;
imagesc(data{imageNumBright})
colormap('gray')
clim([0 13000])
set(gca,'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')
exportgraphics(gca,fullfile(saveDir,'figS6a_brightfield_raw.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Raw fluorescence image
figure;
imagesc(data{imageNum})
colormap('gray')
clim([1000 12000])
set(gca,'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')
exportgraphics(gca,fullfile(saveDir,'figS6a_fluorescence_raw.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Raw fluoresence surface plot
hold on
surf(data{imageNum},'Linestyle','none')
clim([1000 17000])
set(gca,'XLim',[1 size(data{1},2)],'YLim',[1 size(data{1},1)],'ZLim',[0 13000],...
    'View',[-26 8],'Box','on','XTick',[],'YTick',[],'ZTick',[])
exportgraphics(gca,fullfile(saveDir,'figS6a_fluorescence_raw_surf.png'),...
    'ContentType','image','Resolution',300);

%% Figure S6B - Example background image

indBackground = backgroundInfo.imageNum == imageNumBack;

% Raw background image
figure
imagesc(data{imageNumBack})
colormap('gray')
clim([1000 6000])
set(gca,'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')
exportgraphics(gca,fullfile(saveDir,'figS6b_background_raw.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Raw background surface plot
hold on
surf(data{imageNumBack},'Linestyle','none')
clim([1000 10000])
set(gca,'XLim',[1 size(data{1},2)],'YLim',[1 size(data{1},1)],'ZLim',[0 13000],...
    'View',[-26 8],'Box','on','XTick',[],'YTick',[],'ZTick',[])
exportgraphics(gca,fullfile(saveDir,'figS6b_background_raw_surf.png'),...
    'ContentType','image','Resolution',300);

% Smoothed background surface plot
figure
imagesc(background{indBackground})
colormap('gray')
clim([1000 6000])
set(gca,'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')
hold on
surf(background{indBackground},'Linestyle','none')
clim([1000 10000])
set(gca,'XLim',[1 size(data{1},2)],'YLim',[1 size(data{1},1)],'ZLim',[0 13000],...
    'View',[-26 8],'Box','on','XTick',[],'YTick',[],'ZTick',[])
exportgraphics(gca,fullfile(saveDir,'figS6b_background_smooth_surf.png'),...
    'ContentType','image','Resolution',300);

%% Figure S6C - Normalized Fluorescent Images

indProfile = lawnProfiles.imageNum == imageNum;

% Normalized fluorescence image
figure;
imagesc(lawnProfiles.imageNormalized{indProfile})
colormap('gray')
clim([1000 12000])
set(gca,'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')
exportgraphics(gca,fullfile(saveDir,'figS6c_fluorescence_norm.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Normalized fluorescence surface plot
hold on
surf(lawnProfiles.imageNormalized{indProfile},...
    'Linestyle','none')
clim([1000 17000])
set(gca,'XLim',[1 size(data{1},2)],'YLim',[1 size(data{1},1)],'ZLim',[0 13000],...
    'View',[-26 8],'Box','on','XTick',[],'YTick',[],'ZTick',[])
exportgraphics(gca,fullfile(saveDir,'figS6c_fluorescence_norm_surf.png'),...
    'ContentType','image','Resolution',300);

%% Figure S6D - Detection of patches and extraction of cross-sectional profile

lawnMask = lawnProfiles.labeled{indProfile};

% Lawn and template mask image
figure;
imagesc((bwdist(~lawnMask) > 20 | ~lawnMask) & ...
    (bwdist(~templateMask)>20 | ~templateMask))
for numROI = 1:max(lawnMask,[],'all')
    [x,y] = find(lawnMask == numROI);
    text(mean(y),mean(x),num2str(numROI),'FontSize',12,'FontName','Arial',...
        'VerticalAlignment','middle','HorizontalAlignment','center')
end
colormap('gray')
hold on
plot(size(lawnMask,2)-225-[0 5*scale],ones(2,1).*(size(lawnMask,1) - 240),'k','LineWidth',2)
text(size(lawnMask,2)-225-2.5*scale,size(lawnMask,1)-220,'5 mm','FontSize',12,...
    'FontName','Arial','VerticalAlignment','top','HorizontalAlignment','center')
set(gca,'XTick',[],'YTick',[])
set(gcf,'Toolbar','none')
exportgraphics(gca,fullfile(saveDir,'figS6d_lawn_template_mask.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% ROI image
figure('Position',[840 638 420 420]);
ROImask = lawnProfiles.imageNormalized{indProfile}.*(lawnMask == 1);
[x,y] = find(ROImask);
ROImask(ROImask == 0) = 12000;
imagesc(ROImask((mean(x)-130):(mean(x)+130),(mean(y)-130):(mean(y)+130)))
colormap('gray')
clim([1000 12000])
hold on
scatter(130.5,130.5,'k')
plot(255-[0 scale/2],ones(2,1).*235,'k','LineWidth',1)
text(255-scale/4,240,'500 μm','FontSize',8,...
    'FontName','Arial','VerticalAlignment','top','HorizontalAlignment','center')
set(gca,'XTick',[],'YTick',[],'Units','inches')
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 1.5])
exportgraphics(gca,fullfile(saveDir,'figS6d_ROI_mask.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% ROI lawn Profile
figure;
plot(lawnProfiles.distances{indProfile},...
    lawnProfiles.pixelValuesNormalized{indProfile}(1,:),'k','LineWidth',1.5)
set(gca,'FontName','Arial','FontSize',8,'XLim',[-0.55 0.85],'XTick',-0.5:0.25:0.75,'Units','inches')
hold on
xlabel('Distance from patch border (μm)','FontSize',8)
ylabel('Pixel intensity (a.u.)','FontSize',8)
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 1.25])
exportgraphics(gca,fullfile(saveDir,'figS6d_ROI_profile.pdf'),...
    'BackgroundColor','none','ContentType','vector')

%% Figure S6E - Border Amplitude Detection + Linear Fit over Time

indLawn = find(lawnAnalysis.imageNum == imageNum,1,'first');

% ROI lawn Profile
figure; hold on
plot(lawnProfiles.distances{indProfile},...
    lawnProfiles.pixelValuesNormalized{indProfile}(1,:),...
    'k','LineWidth',1.5)
scatter(lawnAnalysis.xOuterEdge(indLawn),lawnAnalysis.yOuterEdge(indLawn),'b')
scatter(lawnAnalysis.xPeak(indLawn),lawnAnalysis.yPeak(indLawn),'b')
plot([-0.1 -0.1],[lawnAnalysis.yOuterEdge(indLawn) lawnAnalysis.yPeak(indLawn)],'b','LineWidth',1)
plot([-0.15 -0.05],[lawnAnalysis.yOuterEdge(indLawn) lawnAnalysis.yOuterEdge(indLawn)],'b','LineWidth',1)
plot([-0.15 -0.05],[lawnAnalysis.yPeak(indLawn) lawnAnalysis.yPeak(indLawn)],'b','LineWidth',1)
text(-0.3,mean([lawnAnalysis.yOuterEdge(indLawn) lawnAnalysis.yPeak(indLawn)]),sprintf('border\namplitude'),...
    'Color','b','FontSize',8,'FontName','Arial','VerticalAlignment','middle','HorizontalAlignment','center')
indLawn = find(lawnAnalysis.imageNum == imageNum,1,'first');
set(gca,'FontName','Arial','FontSize',8,'XLim',[-0.55 0.85],'XTick',-0.5:0.25:0.75,'Units','inches','Box','on')
hold on
xlabel('Distance from patch border (μm)','FontSize',8)
ylabel('Pixel intensity (a.u.)','FontSize',8)
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 1.25])
exportgraphics(gca,fullfile(saveDir,'figS6e_ROI_profile.pdf'),...
    'BackgroundColor','none','ContentType','vector')
%%
% Example Border Amplitude Over Time
conditionNum = find(strcmp(conditions.peptone,'with') & conditions.lawnVolume == 0.5 & ...
    conditions.growthTimeCondition == 0 & conditions.OD600 == 10);
indAnalysis = lawnAnalysis.condition == conditionNum & lawnAnalysis.include;

tFit = (0:360) + 60*borderAmplitude.growthTimeCondition;
indBorderAmp = find(strcmp(borderAmplitude.peptone,'with') & borderAmplitude.lawnVolume == 0.5 & borderAmplitude.OD600 == 10);

figure; hold on

gscatter(minutes(lawnAnalysis.growthTimeTotal(indAnalysis)),...
    lawnAnalysis.borderAmplitude(indAnalysis)./lawnAnalysis.exposureTime(indAnalysis),...
    lawnAnalysis.expNum(indAnalysis))
plot(tFit(indBorderAmp,:)-60*borderAmplitude.growthTimeCondition(indBorderAmp),...
    tFit(indBorderAmp,:).*borderAmplitude.slope(indBorderAmp) + borderAmplitude.intercept(indBorderAmp),...
    'k','LineWidth',1)
legend('off')
set(gca,'FontName','Arial','FontSize',8,'XLim',[0 240],'XTick',0:60:240,'YLim',[0 5],...
    'Units','inches','Box','on')
xlabel('time (min)','FontSize',8)
ylabel('relative border amplitude (a.u.)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 2 1.25])
exportgraphics(gca,fullfile(saveDir,'figS6e_borderAmplitude.pdf'),...
    'BackgroundColor','none','ContentType','vector')