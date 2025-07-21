%% Load data into workspace

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'encounter.mat'),'encounter');
% borderAmplitude = readtable([path,'foragingGFP\borderAmplitude.csv']);

saveDir = [path,'figures\FigS15\'];

rhoKColor = [61 107 102]./255; % turquoise
tauSColor = [178 73 36]./255; % rust
rhoHColor = [107 61 90]./255; % marroon
rhoEColor = [69 61 109]./255; % indigo
permuteColor = [204 121 167]./255; % pink

% Get conditions
[G,GID] = findgroups(encounter(:,{'expName','lawnVolume',...
    'growthCondition','OD600Label','strainName','strainID'}));
conditionG = find(strcmp(GID.expName,expName) & GID.lawnVolume == 0.5 & ...
    ~(strcmp(GID.OD600Label,'0.10') & GID.growthCondition == 48));

%% Figure S15A - Plot Log-Likelihood, AIC, BIC of all GLM combinations

load(fullfile(path,'geometricData','geometricData_combinations.mat'));
nModels = length(beta);

figure('Position',[100 100 1000 800]);
for i = 1:3
    subplot(414);
    for m = 1:nModels
        switch i
            case 1
                v = Violin({logLike(m,:,:)-logLike(end,:,:)},m,'ViolinColor',{[0 0 0]},...
                    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
                ylim([-300 25]); ylabel('Δ log-likelihood','FontSize',8)
                [~,bestModel] = max(logLike);
                textY = min(logLike(m,:,:)-logLike(end,:,:),[],'all')-50;
            case 2
                v = Violin({AIC(m,:,:)-AIC(end,:,:)},m,'ViolinColor',{[0 0 0]},...
                    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
                ylim([-50 600]); ylabel('Δ AIC','FontSize',8)
                [~,bestModel] = min(AIC);
                textY = max(AIC(m,:,:)-AIC(end,:,:),[],'all')+100;
            case 3
                v = Violin({BIC(m,:,:)-BIC(end,:,:)},m,'ViolinColor',{[0 0 0]},...
                    'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
                ylim([-50 600]); ylabel('Δ BIC','FontSize',8)
                [~,bestModel] = min(BIC);
                textY = max(BIC(m,:,:)-BIC(end,:,:),[],'all')+100;
        end
        v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
        v.MeanPlot.LineWidth = 0.75;
    end
    yline(0,'k')
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[]);
    set(gca,'Position',[0.9 8-1.15*i 7 1]);
    xlim(0.5 + [0 nModels]);

    bestModels = unique(bestModel);
    bestCount = arrayfun(@(m) sum(bestModel == m,'all'),bestModels);
    for m = 1:length(bestModels)
        text(bestModels(m),textY,[num2str(100*bestCount(m)/sum(bestCount),'%.1f'),'%'],...
            'HorizontalAlignment','center','FontName','Arial','FontSize',8)
    end
end
xticks(1:nModels); xticklabels({''})

exportgraphics(gcf,[saveDir,'figS15a.pdf'],'ContentType','vector')

%% Figure S15B - Plot Log-Likelihood, AIC, BIC of nested GLMs

load(fullfile(path,'geometricData','geometricData.mat'));
allVars = vertcat(modelVars{:});

nModels = length(beta);
nReps = length(modelVars);
nBoot = size(bootWorm,2);
betaNames = {'β_0','+ β_k ρ_k','+ β_s τ_s','+ β_h ρ_h','+ β_e ρ_e'};

figure('Position',[100 100 800 600]);
colormap(jet)
for m = 1:nModels
    subplot(3,nModels,3*nModels); imagesc(squeeze(logLike(m,1:nReps,1:nBoot)))
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 3.8 0.75 0.75]); clim(prctile(logLike(:),[0 100]));
    title(betaNames{m},'FontSize',9,'FontWeight','normal','FontName','Times New Roman','FontAngle','italic');

    if m == nModels
        c = colorbar('Units','inches','Position',[(nModels+1)*0.9 3.8 0.1 0.75]);
        clim(prctile(logLike(:),[0 100])); c.Label.String = 'log-likelihood'; c.Label.Rotation = -90;
    end

    subplot(3,nModels,3*nModels); imagesc(squeeze(AIC(m,1:nReps,1:nBoot)))
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 2.9 0.75 0.75]); clim(prctile(AIC(:),[0 100]));

    if m == nModels
        c = colorbar('Units','inches','Position',[(nModels+1)*0.9 2.9 0.1 0.75]);
        clim(prctile(AIC(:),[0 100])); c.Label.String = 'AIC'; c.Label.Rotation = -90;
    end
    
    subplot(3,nModels,3*nModels); imagesc(squeeze(BIC(m,1:nReps,1:nBoot)))
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 2 0.75 0.75]); clim(prctile(AIC(:),[0 100]));

    if m == 1
        set(gca,'XTick',[1 nBoot],'YTick',[1 nReps])
        xlabel({'worm','samples'},'FontSize',8','FontName','Arial');
        ylabel({'encounter','samples'},'FontSize',8','FontName','Arial');
    elseif m == nModels
        c = colorbar('Units','inches','Position',[(nModels+1)*0.9 2 0.1 0.75]);
        clim(prctile(BIC(:),[0 100])); c.Label.String = 'BIC'; c.Label.Rotation = -90;
    end
end

exportgraphics(gcf,[saveDir,'figS15b.pdf'],'ContentType','vector')

%% Figure S15C - Plot beta values for shuffles

betaColors = [0 0 0; rhoKColor; tauSColor; rhoHColor; rhoEColor];

nModels = length(beta);
bestModel = nModels;
meanBeta = mean(betaShuffle,[2,3]);
pBeta = mean(betaShuffle.*sign(meanBeta) < 0,[2,3])
% betaCI = arrayfun(@(n) coefCI(mdl{bestModel,n}),1:nReps);
% betaShuffleCI = prctile(betaShuffle,[2.5 50 97.5],2);

figure; hold on
yline(0,'k')
for m = 1:nModels
    % v = Violin({betaShuffle(m,:,:)},m,'ViolinColor',{permuteColor},...
    %     'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
    % v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;

    v = Violin({betaShuffle(m,:,:)},m,'ViolinColor',{[0 0 0]},...
        'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
    v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
    v.EdgeColor = betaColors(m,:); v.WhiskerPlot.Color = betaColors(m,:);
    v.MedianPlot.MarkerEdgeColor = betaColors(m,:); v.MeanPlot.Color = betaColors(m,:);
    v.MeanPlot.LineWidth = 0.75;
    
    % [~,pBeta(m)] = ttest2(reshape(beta{bestModel}(m,:,:),[],1),...
    %     reshape(betaShuffle(m,:,:),[],1));

    if pBeta(m) < 0.001
            sig = '***';
        elseif pBeta(m) < 0.01
            sig = '**';
        elseif pBeta(m) < 0.05
            sig = '*';
        else
            sig = 'n.s.';
        end
    text(m,4.5,sig,'Color',betaColors(m,:),'VerticalAlignment','bottom',...
        'HorizontalAlignment','center','FontSize',8,'FontName','Arial')
end
ylim([-2.25 5.25]); xlim(0.5 + [0 nModels]); xticks(1:nModels); xticklabels({'1','ρ_k','τ_s','ρ_h','ρ_e'})
set(gca,'Units','inches','FontName','Arial','FontSize',8,...
    'Box','on','TickDir','out');
ylabel('β^*','FontSize',9,'FontName','Times New Roman','FontAngle','italic')
ax = gca; ax.XAxis.FontName = 'Times New Roman'; ax.XAxis.FontAngle = 'italic'; ax.XAxis.FontSize = 9;
set(gca,'Position',[ax.Position(1:2) 0.975 1.75])
exportgraphics(gcf,[saveDir,'figS15c.pdf'],'ContentType','vector')

%% Plot beta values across repetitions

betaNames = {'β_0','β_k','β_s','β_h','β_e'};
col = [linspace(0,1,100)',linspace(0,1,100)',ones(100,1);...
    ones(100,1),linspace(1,0,100)',linspace(1,0,100)'];
betaMaxMin = prctile(beta{end},[0 100],'all');
betaLim = max(abs(betaMaxMin))*[-1 1];

figure('Position',[100 100 800 600]);
colormap(col)
for m = 1:nModels    
    subplot(2,nModels,nModels); imagesc(squeeze(beta{end}(m,:,:)))
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 2.9 0.75 0.75]); clim(betaLim);
    title(betaNames{m},'FontSize',9,'FontWeight','bold','FontName',...
        'Times New Roman','FontAngle','italic','Color',betaColors(m,:));

     if m == nModels
        c = colorbar('Units','inches','Position',[(nModels+1)*0.9 2.9 0.1 0.75]);
        clim(betaLim); c.Label.String = 'β'; c.Label.Rotation = 0;
        c.Label.FontName = 'Times New Roman'; c.Label.FontAngle = 'italic';
        c.Label.FontWeight = 'bold';
    end

    subplot(2,nModels,nModels); imagesc(squeeze(betaShuffle(m,:,:)))
    set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out','XTick',[],'YTick',[]);
    set(gca,'Position',[m*0.9 2 0.75 0.75]); clim(betaLim);

    if m == 1
        set(gca,'XTick',[1 nBoot],'YTick',[1 nReps])
        xlabel({'worm','samples'},'FontSize',8','FontName','Arial');
        ylabel({'encounter','samples'},'FontSize',8','FontName','Arial');
    elseif m == nModels
        c = colorbar('Units','inches','Position',[(nModels+1)*0.9 2 0.1 0.75]);
        clim(betaLim); c.Label.String = 'β_{shuffled}'; c.Label.Rotation = 0;
        c.Label.FontName = 'Times New Roman'; c.Label.FontAngle = 'italic';
        c.Label.FontWeight = 'bold'; c.Label.Color = permuteColor
    end
    
end

% exportgraphics(gcf,[saveDir,'figS15b.pdf'],'ContentType','vector')
