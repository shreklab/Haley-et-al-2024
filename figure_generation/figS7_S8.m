% Load data

path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'foragingGFP','analyzeGFP_24-04-04.mat'),...
    'background','backgroundInfo','data','info',...
    'lawnAnalysis','lawnProfiles','metaData','template');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
saveDir = [path,'figures\FigS7_S8\'];

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

% Get conditions
[lawnAnalysis.condition,conditions] = ...
    findgroups(lawnAnalysis(:,{'peptone','lawnVolume','growthTimeCondition',...
    'OD600'}));
condToAnalyze = find(ismember(conditions.growthTimeCondition,[0,48]) | ...
    (conditions.growthTimeCondition == 12 & conditions.OD600 == 1 & conditions.lawnVolume == 0.5));

% Get color values for plotting
borderValue = borderAmplitude.slope.*60.*(borderAmplitude.growthTimeCondition+1) + ...
        borderAmplitude.intercept;
borderAmp10 = borderValue(strcmp(borderAmplitude.peptone,'with') & ...
    borderAmplitude.OD600 == 10.00 & borderAmplitude.lawnVolume == 0.5);
relativeBorder = 10.*borderValue./borderAmp10;
colorValue = min(1 - (log(relativeBorder)*0.09 + 0.4),0.8);

%% Figure S7A - Patch examples for most concentrated peptone+ conditions

% expNum = 3+5, expTime = 2 hr
imageNums = [749,649,660,2141,2334,2527,2726,2903,3077,1147,521]; % 10,12H,48H,1(20),2(20),3(20),4(20),5(20),10(20),48H(20),48H(200)
patchNums = [1,2,1,1,1,1,1,1,1,1,1];
lineColor = repmat(linspace(0.8,0,12)',1,3); lineColor = lineColor([10:12,12],:);
lineWidth = 0.3*([10:12,12]);

figure(1)
set(gcf,'Position',[840 638 420 420]);
for i = 1:length(imageNums)
    imageNum = imageNums(i);
    lawnMask = lawnProfiles.labeled{lawnProfiles.imageNum == imageNum};
    [x,y] = find(lawnMask == patchNums(i));
    ROImask = lawnProfiles.imageNormalized{lawnProfiles.imageNum == imageNum};
    indLawn = find(lawnAnalysis.imageNum == imageNum);
    indLawn = indLawn(patchNums(i));
    if lawnAnalysis.lawnVolume(indLawn) == 0.5
        offset = 130;
        xMax = 1.2;
    else
        offset = 130*8;
        xMax = 10.2;
        figure(1)
        set(gcf,'Position',[0 0 420*2 420*2])
        figure(2)
        set(gcf,'Position',[0 0 560*3 420])
    end

    figure(1)
    imagesc(ROImask((mean(x)-offset):(mean(x)+offset),...
        (mean(y)-offset):(mean(y)+offset))./lawnAnalysis.exposureTime(indLawn))
    colormap(flipud(gray))
    clim(lawnAnalysis.yOuterEdge(indLawn)/lawnAnalysis.exposureTime(indLawn) + ...
        [0 49.2])
    set(gca,'XTick',[],'YTick',[],'Position',[0.13,0.11,0.775,0.775])
    set(gcf,'Toolbar','none')
    exportgraphics(gca,fullfile(saveDir,...
        ['figS7a_ROI_OD',num2str(100*lawnAnalysis.OD600(indLawn),'%.4d'),'_',...
        num2str(lawnAnalysis.growthTimeCondition(indLawn),'%.2d'),'H_',...
        num2str(10*lawnAnalysis.lawnVolume(indLawn),'%.4d'),'uL_high.pdf']),...
        'BackgroundColor','none','ContentType','vector')

    % if lawnAnalysis.lawnVolume(indLawn) == 0.5
    %     figure(3); hold on
    % elseif lawnAnalysis.lawnVolume(indLawn) == 20
    %     figure(4); hold on
    % else
    %     figure(5); hold on
    % end
    % indColor = strcmp(borderAmplitude.peptone,lawnAnalysis.peptone{indLawn}) & ...
    %     borderAmplitude.lawnVolume == lawnAnalysis.lawnVolume(indLawn) & ...
    %     borderAmplitude.OD600 == lawnAnalysis.OD600(indLawn) & ...
    %     borderAmplitude.growthTimeCondition == lawnAnalysis.growthTimeCondition(indLawn);
    % x = lawnProfiles.distances{lawnProfiles.imageNum == imageNum};
    % y = (lawnProfiles.pixelValuesNormalized{lawnProfiles.imageNum == imageNum} - ...
    %     lawnAnalysis.yOuterEdge(indLawn))./lawnAnalysis.exposureTime(indLawn);
    % plot(x,y(patchNums(i),:),'Color',colorValue(indColor).*[1 1 1],'LineWidth',abs(1-colorValue(indColor))*4)
end

% figure(3)
% set(gca,'FontName','Arial','FontSize',8,'Units','inches',...
%     'XLim',[-0.5 1.2],'XTick',0:1,'YLim',[-5 70],'Box','on')
% ax = gca; set(gca,'Position',[ax.Position(1:2) 4*(1.7/10.5) 0.75])
% xlabel('Distance from patch border (μm)','FontSize',8)
% ylabel('Relative Fluorescence','FontSize',8)
% exportgraphics(gca,fullfile(saveDir,'figS7a_profile_0005uL_high.pdf'),...
%     'BackgroundColor','none','ContentType','vector')
% 
% figure(4)
% set(gca,'FontName','Arial','FontSize',8,'Units','inches',...
%     'XLim',[-0.5 4.5],'XTick',0:4,'YLim',[-5 70],'Box','on')
% ax = gca; set(gca,'Position',[ax.Position(1:2) 4*(5/10.5) 0.75])
% xlabel('Distance from patch border (μm)','FontSize',8)
% ylabel('Relative Fluorescence','FontSize',8)
% exportgraphics(gca,fullfile(saveDir,'figS7a_profile_0200uL_high.pdf'),...
%     'BackgroundColor','none','ContentType','vector')
% 
% figure(5)
% set(gca,'FontName','Arial','FontSize',8,'Units','inches',...
%     'XLim',[-0.5 10],'XTick',0:10,'YLim',[-5 70],'Box','on')
% ax = gca; set(gca,'Position',[ax.Position(1:2) 4 0.75])
% xlabel('Distance from patch border (μm)','FontSize',8)
% ylabel('Relative Fluorescence','FontSize',8)
% exportgraphics(gca,fullfile(saveDir,'figS7a_profile_2000uL_high.pdf'),...
%     'BackgroundColor','none','ContentType','vector')


%% Figure S7A - Patch examples for least concentrated peptone+ conditions

% expNum = 3, expTime = 2 hr
imageNums = [615,544,575,584,535,595,555]; % 0.5,1,2,3,4,5,10
patchNums = [1,2,3,1,1,2,3];

figure(1)
set(gcf,'Position',[840 638 420 420]);

for i = 1:length(imageNums)
    imageNum = imageNums(i);
    lawnMask = lawnProfiles.labeled{lawnProfiles.imageNum == imageNum};
    [x,y] = find(lawnMask == patchNums(i));
    ROImask = lawnProfiles.imageNormalized{lawnProfiles.imageNum == imageNum};
    indLawn = find(lawnAnalysis.imageNum == imageNum);
    indLawn = indLawn(patchNums(i));
    if lawnAnalysis.lawnVolume(indLawn) == 0.5
        offset = 130;
        xMax = 1.2;
    else
        offset = 130*8;
        xMax = 10.2;
        figure(1)
        set(gcf,'Position',[0 0 420*2 420*2])
        figure(2)
        set(gcf,'Position',[0 0 560*3 420])
    end

    figure(1)
    imagesc(ROImask((mean(x)-offset):(mean(x)+offset),...
        (mean(y)-offset):(mean(y)+offset))./lawnAnalysis.exposureTime(indLawn))
    colormap(flipud(gray))
    clim(lawnAnalysis.yOuterEdge(indLawn)/lawnAnalysis.exposureTime(indLawn) + ...
        [0 2.2])
    set(gca,'XTick',[],'YTick',[],'Position',[0.13,0.11,0.775,0.775])
    set(gcf,'Toolbar','none')
    exportgraphics(gca,fullfile(saveDir,...
        ['figS7a_ROI_OD',num2str(100*lawnAnalysis.OD600(indLawn),'%.4d'),'_',...
        num2str(lawnAnalysis.growthTimeCondition(indLawn),'%.2d'),'H_',...
        num2str(10*lawnAnalysis.lawnVolume(indLawn),'%.4d'),'uL_noPeptone.pdf']),...
        'BackgroundColor','none','ContentType','vector')
end

%% Figure S7A - Patch examples for peptone- conditions

% expNum = 3, expTime = 2 hr
imageNums = [740,789,700,780,760,729,749,1925,2141]; % 0.5,1,2,3,4,5,10,0.5(20),1(20)
patchNums = [2,1,1,1,1,11,1,1,1];

figure(1)
set(gcf,'Position',[840 638 420 420]);

for i = 1:length(imageNums)
    imageNum = imageNums(i);
    lawnMask = lawnProfiles.labeled{lawnProfiles.imageNum == imageNum};
    [x,y] = find(lawnMask == patchNums(i));
    ROImask = lawnProfiles.imageNormalized{lawnProfiles.imageNum == imageNum};
    indLawn = find(lawnAnalysis.imageNum == imageNum);
    indLawn = indLawn(patchNums(i));
    if lawnAnalysis.lawnVolume(indLawn) == 0.5
        offset = 130;
        xMax = 1.2;
    else
        offset = 130*8;
        xMax = 10.2;
        figure(1)
        set(gcf,'Position',[0 0 420*2 420*2])
        figure(2)
        set(gcf,'Position',[0 0 560*3 420])
    end

    figure(1)
    imagesc(ROImask((mean(x)-offset):(mean(x)+offset),...
        (mean(y)-offset):(mean(y)+offset))./lawnAnalysis.exposureTime(indLawn))
    colormap(flipud(gray))
    clim(lawnAnalysis.yOuterEdge(indLawn)/lawnAnalysis.exposureTime(indLawn) + ...
        [0 2.2])
    set(gca,'XTick',[],'YTick',[],'Position',[0.13,0.11,0.775,0.775])
    set(gcf,'Toolbar','none')
    exportgraphics(gca,fullfile(saveDir,...
        ['figS7a_ROI_OD',num2str(100*lawnAnalysis.OD600(indLawn),'%.4d'),'_',...
        num2str(lawnAnalysis.growthTimeCondition(indLawn),'%.2d'),'H_',...
        num2str(10*lawnAnalysis.lawnVolume(indLawn),'%.4d'),'uL.pdf']),...
        'BackgroundColor','none','ContentType','vector')
end

%%
conditions.conditionNum(:) =  1:height(conditions);
conditionsPlot = conditions(condToAnalyze,:);
conditionsPlot = sortrows(conditionsPlot,{'peptone','lawnVolume','growthTimeCondition','OD600',},{'ascend','descend','ascend','ascend'});
xOffset = cumsum([0;any([diff(conditionsPlot.lawnVolume),diff(strcmp(conditionsPlot.peptone,'with'))] ~= 0,2)]);

violinHalf = repmat({'full'},length(condToAnalyze),1);
metric = 'borderAmplitude';

figure('Position',[360 278 840 420]);
hold on
for j = 1:length(condToAnalyze)
    ind = lawnAnalysis.include & lawnAnalysis.condition == conditionsPlot.conditionNum(j) & ...
        hours(lawnAnalysis.growthTimeTotal) <= 4 + conditionsPlot.growthTimeCondition(j);
    metricData = 10.*lawnAnalysis.(metric)(ind)./(lawnAnalysis.exposureTime(ind)*borderAmp10);
    indColor = strcmp(borderAmplitude.peptone,conditionsPlot.peptone{j}) & ...
        borderAmplitude.lawnVolume == conditionsPlot.lawnVolume(j) & ...
        borderAmplitude.OD600 == conditionsPlot.OD600(j) & ...
        borderAmplitude.growthTimeCondition == conditionsPlot.growthTimeCondition(j);
    Violin({metricData},j+xOffset(j),'ViolinColor',...
        {colorValue(indColor).*[1 1 1]},'QuartileStyle','shadow','ShowMean',true,...
        'HalfViolin',violinHalf{j});
end
xlim([0 length(condToAnalyze)+1+max(xOffset)])
xticks((1:length(condToAnalyze))+xOffset'); xticklabels(conditionsPlot.OD600)
ylabel('relative density','FontSize',8); ylim([0.05 1600]);
yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'TickDir','out','FontName','Arial','FontSize',8,'LineWidth',1,...
    'YScale','log','Units','inches','Box','on')
set(gcf,'Toolbar','none')
ax = gca; set(gca,'Position',[ax.Position(1:2) 7 1.5])
exportgraphics(gca,fullfile(saveDir,'figS7b_borderAmplitude.pdf'),...
    'BackgroundColor','none','ContentType','vector')

%% Plot fits of border amplitude over time

% Plot 0.5 uL patches w/o peptone
figure('Position',[400 400 1000 300]);
ind = find(strcmp(borderAmplitude.peptone,'without') & borderAmplitude.lawnVolume == 0.5);
for i = 1:length(ind)
    hold on
    plot(tFit(ind(i),:)-60*borderAmplitude.growthTimeCondition(ind(i)),...
        10*(tFit(ind(i),:).*borderAmplitude.slope(ind(i)) + borderAmplitude.intercept(ind(i)))./borderAmp10,...
        'Color',colorValue(ind(i),:).*[1 1 1],'LineWidth',2)
end
xlim([0 360]); xticks(0:60:360)
ylim([0.05 1600]); yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'YScale','log','Box','on','FontName','Arial','FontSize',8,'LineWidth',1,'Units','inches','TickDir','Out')
xlabel('experiment time (min)','FontSize',8)
ylabel('relative density','FontSize',8)
title('0.5 μL (no peptone)','FontSize',8,'FontWeight','normal','FontName','Arial')
set(gca,'Position',[1 1 2 1.5])
exportgraphics(gca,fullfile(saveDir,'figS7c_borderAmplitude_0005uL_noPeptone.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Plot 0.5 uL patches w/ peptone
figure('Position',[400 400 1000 300]);
ind = find(strcmp(borderAmplitude.peptone,'with') & borderAmplitude.lawnVolume == 0.5);
for i = 1:length(ind)
    hold on
    plot(tFit(ind(i),:)-60*borderAmplitude.growthTimeCondition(ind(i)),...
        10*(tFit(ind(i),:).*borderAmplitude.slope(ind(i)) + borderAmplitude.intercept(ind(i)))./borderAmp10,...
        'Color',colorValue(ind(i),:).*[1 1 1],'LineWidth',2)
end
xlim([0 360]); xticks(0:60:360)
ylim([0.05 1600]); yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'YScale','log','Box','on','FontName','Arial','FontSize',8,'LineWidth',1,'Units','inches','TickDir','Out')
xlabel('experiment time (min)','FontSize',8)
ylabel('relative density','FontSize',8)
title('0.5 μL','FontSize',8,'FontWeight','normal','FontName','Arial')
set(gca,'Position',[1 1 2 1.5])
exportgraphics(gca,fullfile(saveDir,'figS7c_borderAmplitude_0005uL.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Plot 20 uL patches w/ peptone
figure('Position',[400 400 1000 300]);
ind = find(strcmp(borderAmplitude.peptone,'with') & borderAmplitude.lawnVolume == 20);
for i = 1:length(ind)
    hold on
    plot(tFit(ind(i),:)-60*borderAmplitude.growthTimeCondition(ind(i)),...
        10*(tFit(ind(i),:).*borderAmplitude.slope(ind(i)) + borderAmplitude.intercept(ind(i)))./borderAmp10,...
        'Color',colorValue(ind(i),:).*[1 1 1],'LineWidth',2)
end
xlim([0 360]); xticks(0:60:360)
ylim([0.05 1600]); yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'YScale','log','Box','on','FontName','Arial','FontSize',8,'LineWidth',1,'Units','inches','TickDir','Out')
xlabel('experiment time (min)','FontSize',8)
ylabel('relative density','FontSize',8)
title('20 μL','FontSize',8,'FontWeight','normal','FontName','Arial')
set(gca,'Position',[1 1 2 1.5])
exportgraphics(gca,fullfile(saveDir,'figS7c_borderAmplitude_0200uL.pdf'),...
    'BackgroundColor','none','ContentType','vector')

% Plot 200 uL patches w/ peptone
figure('Position',[400 400 1000 300]);
ind = find(strcmp(borderAmplitude.peptone,'with') & borderAmplitude.lawnVolume == 200);
for i = 1:length(ind)
    hold on
    plot(tFit(ind(i),:)-60*borderAmplitude.growthTimeCondition(ind(i)),...
        10*(tFit(ind(i),:).*borderAmplitude.slope(ind(i)) + borderAmplitude.intercept(ind(i)))./borderAmp10,...
        'Color',colorValue(ind(i),:).*[1 1 1],'LineWidth',2)
end
xlim([0 360]); xticks(0:60:360)
ylim([0.05 1600]); yticks([0.1 1 10 100 1000]); yticklabels([0.1 1 10 100 1000])
set(gca,'YScale','log','Box','on','FontName','Arial','FontSize',8,'LineWidth',1,'Units','inches','TickDir','Out')
xlabel('experiment time (min)','FontSize',8)
ylabel('relative density','FontSize',8)
title('200 μL','FontSize',8,'FontWeight','normal','FontName','Arial')
set(gca,'Position',[1 1 2 1.5])
exportgraphics(gca,fullfile(saveDir,'figS7c_borderAmplitude_2000uL.pdf'),...
    'BackgroundColor','none','ContentType','vector')
