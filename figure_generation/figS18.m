%% Load data into workspace

expName = 'foragingMatching';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,'geometricData','geometricData.mat'));
load(fullfile(path,'encounter.mat'),'encounter');
bodyPart = 'midpoint';
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);
saveDir = [path,'figures\FigS18\'];

exploitColor = [0 114 178]./255; % blue
sampleColor = [0 158 115]./255; % green
searchOnColor = [230 159 0]./255; % orange
rhoKColor = [61 107 102]./255; % turquoise
tauSColor = [178 73 36]./255; % rust
betaColors = [0.5 0.5 0.5; rhoKColor; tauSColor];

%% Get conditions to analyze

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));

conditionG = find(strcmp(GID.expName,expName));
[~,cOrder] = sort(cellfun(@str2num,erase(erase(GID.OD600Label(conditionG),'.'),' ')));
conditionG = conditionG(cOrder);
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
wormGroupInd = arrayfun(@(g) wormGroup == g,conditionG,'UniformOutput',false)';
[~,wormGroupInd] = max([wormGroupInd{:}],[],2); % continuous group id
% [~,~,wormGroupInd] = unique(wormGroup); % continuous group id

% Get indices of encounters labeled exploit, sample, or searchOn
indEncounter = strcmp(encounter.expName,expName) & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & ~strcmp(encounter.label,'') & ~strcmp(encounter.label,'searchOff');

% Get relative border amplitude for each worm
wormBorder = arrayfun(@(w) max([10*min(encounter.borderAmplitude(encounter.wormNum == w & ...
    ismember(G,conditionG) & indEncounter))/borderAmp10,borderAmp0]),wormNums);

%%  Figure S18A - Plot example traces

wormNum = [175,166,163,178,191,182,172]; % ['1', '5', '10', '1 5', '1 10', '5 10', '1 5 10']
for i = 1:length(wormNum)
    % plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
end

%% FigS18B - Encounters labeled (with density) for all worms

encounterColors = max(10*encounter.borderAmplitude./borderAmp10,borderAmp0);
encounterColors = min(1 - (log(encounterColors)*0.09 + 0.4),0.8).*ones(1,3);

wormsPerGroup = histcounts(wormGroupInd,(0:numGroups)+0.5)
figure('Position',[500 500 560 800]);
for k = 1:numGroups
    subplot(numGroups,1,1)
    w = wormNums(wormGroupInd == k);
    for i = 1:length(w)
        ind = find(encounter.wormNum == w(i) & strcmp(encounter.expName,expName) & ...
            ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
            strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn')));
        for j = 1:length(ind)
            patch([encounter.timeEnter([ind(j),ind(j)]);encounter.timeExit([ind(j),ind(j)])]./60,...
                i + [-1 0 0 -1]',encounterColors(ind(j),:),'EdgeAlpha',0);
        end
    end
    set(gca,'Units','inches','Visible','on','XTick',[],...
        'FontName','Arial','FontSize',8,'Box','on',...
        'GridLineWidth',0.75,'XLim',[0 120],'YLim',[0 length(w)],...
        'YTick',length(w)/2,'YTickLabel',length(w),'TickDir','in')
    subplotHeight = 0.7*length(w)/max(wormsPerGroup);
    ax = gca; set(gca,'Position',[ax.Position(1), ...
        (numGroups-k+1)*5.5/numGroups + (0.7-subplotHeight)/2, 3, subplotHeight])
end
xlabel('time (min)','FontSize',8)
set(gca,'XTick',0:20:120)
exportgraphics(gcf,[saveDir,'figS18b.pdf'],'ContentType','vector')

%% FigS18C - Encounters labeled (with state) for all worms

encounterColors = encounter.exploitPosterior.*exploitColor + ...
            (1-encounter.exploitPosterior).*encounter.sensePosterior.*sampleColor + ...
            (1-encounter.exploitPosterior).*(1-encounter.sensePosterior).*searchOnColor;
wormsPerGroup = histcounts(wormGroupInd,(0:numGroups)+0.5)
figure('Position',[500 500 560 800]);
for k = 1:numGroups
    subplot(numGroups,1,1)
    w = wormNums(wormGroupInd == k);
    for i = 1:length(w)
        ind = find(encounter.wormNum == w(i) & strcmp(encounter.expName,expName) & ...
            ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
            strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn')));
        for j = 1:length(ind)
            patch([encounter.timeEnter([ind(j),ind(j)]);encounter.timeExit([ind(j),ind(j)])]./60,...
                i + [-1 0 0 -1]',encounterColors(ind(j),:),'EdgeAlpha',0);
        end
    end
    set(gca,'Units','inches','Visible','on','XTick',[],...
        'FontName','Arial','FontSize',8,'Box','on',...
        'GridLineWidth',0.75,'XLim',[0 120],'YLim',[0 length(w)],...
        'YTick',length(w)/2,'YTickLabel',length(w),'TickDir','in')
    subplotHeight = 0.7*length(w)/max(wormsPerGroup);
    ax = gca; set(gca,'Position',[ax.Position(1), ...
        (numGroups-k+1)*5.5/numGroups + (0.7-subplotHeight)/2, 3, subplotHeight])
end
xlabel('time (min)','FontSize',8)
set(gca,'XTick',0:20:120)
exportgraphics(gcf,[saveDir,'figS18c.pdf'],'ContentType','vector')

%% Figure S18E - Plot beta values without history

meanBeta = mean(beta{3},[2,3]);
pBeta = size(beta{3},1)*2*mean(beta{3}.*sign(meanBeta) < 0,[2,3])
indBeta = find(~strcmp(predictors,'timeOffSinceExploit'));

figure; hold on
yline(0,'k')
for m = 1:size(beta{3},1)

    v = Violin({beta{3}(m,:,:)},m,'ViolinColor',{[0 0 0]},...
        'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
    v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
    v.EdgeColor = betaColors(m,:); v.WhiskerPlot.Color = betaColors(m,:);
    v.MedianPlot.MarkerEdgeColor = betaColors(m,:); v.MeanPlot.Color = betaColors(m,:);
    v.MeanPlot.LineWidth = 0.75;

    if pBeta(m) < 0.001
        sig = '***';
    elseif pBeta(m) < 0.01
        sig = '**';
    elseif pBeta(m) < 0.05
        sig = '*';
    else
        sig = 'n.s.';
    end
    text(m,2.25,sig,'Color',betaColors(m,:),'VerticalAlignment','bottom',...
        'HorizontalAlignment','center','FontSize',8,'FontName','Arial')
end
ylim([-3.25 3.25]); xlim(0.5 + [0 size(beta{3},1)]); xticks(1:size(beta{3},1)); xticklabels({'1','ρ_k','τ_s'})
set(gca,'Units','inches','FontName','Arial','FontSize',8,...
    'Box','on','TickDir','out');
ylabel('β','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
set(gca,'Position',[ax.Position(1:2) 0.975 1.25])
exportgraphics(gcf,[saveDir,'figS18e.pdf'],'ContentType','vector')

%% Figure S18D - Beta values of alternative time-dependent term

load(fullfile(path,'geometricData','geometricData_time.mat'));
nModels = length(beta);
combos = {[0:2,4:5];[0:1,3:5];0:5};
betaColors = [0 0 0; rhoKColor; tauSColor; tauSColor; rhoHColor; rhoEColor];
betaLabels = {'1','ρ_k','τ_s','τ_t','ρ_h','ρ_e'};

figure;
for n = 1:nModels
    meanBeta = mean(beta{n},[2,3]);
    pBeta = nModels*2*mean(beta{n}.*sign(meanBeta) < 0,[2,3])
    nTerms = length(pBeta);
    subplot(1,3,n); hold on
    yline(0,'k')
    for m = 1:nTerms

        v = Violin({beta{n}(m,:,:)},m,'ViolinColor',{[0 0 0]},...
            'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
        v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
        v.MeanPlot.LineWidth = 0.75;
        v.EdgeColor = betaColors(combos{n}(m)+1,:); v.WhiskerPlot.Color = betaColors(combos{n}(m)+1,:);
        v.MedianPlot.MarkerEdgeColor = betaColors(combos{n}(m)+1,:); v.MeanPlot.Color = betaColors(combos{n}(m)+1,:);
        v.MeanPlot.LineWidth = 0.75;

        if pBeta(m) < 0.001
            sig = '***';
        elseif pBeta(m) < 0.01
            sig = '**';
        elseif pBeta(m) < 0.05
            sig = '*';
        else
            sig = 'n.s.';
        end
        text(m,4.5,sig,'Color',betaColors(combos{n}(m)+1,:),'VerticalAlignment','bottom',...
            'HorizontalAlignment','center','FontSize',8,'FontName','Arial')
    end
    ylim([-2.25 5.25]); xlim(0.5 + [0 nTerms]); xticks(1:nTerms); xticklabels(betaLabels(combos{n}+1))
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
    if n == 1
        ylabel('β','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
        title('satiety','FontSize',8,'FontWeight','normal','FontName','Arial');
    elseif n == 2
        yticks([])
        title('transfer','FontSize',8,'FontWeight','normal','FontName','Arial');
    else
        yticks([])
        title('satiety and transfer','FontSize',8,'FontWeight','normal','FontName','Arial');
    end
    ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
    set(gca,'Position',[n*1.125 1 0.975*nTerms/5 1.75])
    
end
exportgraphics(gcf,[saveDir,'figS18d.pdf'],'ContentType','vector')