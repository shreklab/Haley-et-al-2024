%% Load Data

path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
load(fullfile(path,'encounter.mat'),'encounter');
load(fullfile(path,'labelEncounters.mat'),'exploitModel','senseModel',...
    'exploitDist','exploreDist','exploitCluster','clusterVar1',...
    'posteriorProbSense','posteriorProbSenseMarginal','minVelocityOnly',...
    'senseCluster','indSense','indNonsense','alpha','posteriorVarExploit','bestModel');
expName = 'foragingConcentration';
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
saveDir = [path,'figures\FigS11_S12\'];

% Define cluster colors
sampleColor = [0 158 115]./255; % green
exploitColor = [0 114 178]./255; % blue
searchColor = [230 159 0]./255; % orange
exploreColor = [213 94 0]./255; % red
colorValue = min(1 - (log(10)*0.09 + 0.4),0.8);

%% Get condition ids

[G,GID] = findgroups(encounter(:,{'expName','lawnVolume','growthCondition',...
    'lawnOD600','peptone'}));
conditionG = find(strcmp(GID.expName,'foragingConcentration') & GID.lawnVolume == 0.5 & ...
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
colorValues = min(1 - (log(relativeBorder)*0.09 + 0.4),0.8);

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
indEncounter = strcmp(encounter.expName,'foragingConcentration') & ismember(encounter.wormNum,wormNums) & ...
    ~encounter.exclude & (strcmp(encounter.label,'exploit') | ...
    strcmp(encounter.label,'sample') | strcmp(encounter.label,'searchOn'));

%% Figure S11A - Heatmap of Avg. Velocity vs. Duration

clusterVars = [log10(encounter.duration./60),log10(encounter.velocityOn)];
include = sum(isnan(clusterVars),2) == 0 & ~encounter.exclude;
clusterVars = clusterVars(include,:);

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
histogram2(clusterVars(:,1),clusterVars(:,2),...
    linspace(-2.1,2.1,51),linspace(1.2,2.8,51),'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','probability');
colormap(flipud(gray))
xlim([-2.1 2.1]); ylim([1.2 2.8]); legend('off')
xticks(log10([2 10 60 600 3600]./60)); xticklabels({'1/30','1/6','1','10','60'});
yticks(log10([25 50 100 200 400])); yticklabels([25 50 100 200 400])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[1 1 3 3])
xlabel('encounter duration (min)','FontSize',8); ylabel('avg. velocity on patch (μm/s)','FontSize',8)

c = colorbar('Units','inches');
set(c,'Position',[4.1 1 0.15 3]);
set(c.Label,'String','p(encounter)','FontSize',8,'Rotation',-90)
exportgraphics(gcf,[saveDir,'figS11a.pdf'],'ContentType','vector')

%% Figure S11B - Projection onto PC1 and Silverman's Critical Bandwidth

% Create mean-centered matrix of variables A
A = clusterVars;
A = A(~any(isnan(A),2),:);
A = A - mean(A,1);

% Compute first eigenvectors of A
[V,~] = eigs(A'*A,1);

% Project values in A onto the eigenvector V
PC1 = A*V(:);

% Get critical bandwidth of data
[criticalBW,silvermansBW,grid] = silvermansTest(PC1,false);

% Plot histogram of values projected onto PC1
figure; hold on;
histogram(PC1,'BinEdges',grid,'Normalization','probability','FaceColor',colorValue.*[1 1 1])
binWidth = grid(2) - grid(1);
p1 = plot(grid,binWidth*mvksdensity(PC1,grid,'Bandwidth',silvermansBW),'LineWidth',2,'Color',exploitColor);
p2 = plot(grid,binWidth*mvksdensity(PC1,grid,'Bandwidth',criticalBW),'LineWidth',2,'Color',exploreColor);
legend([p1,p2],{['Silverman''s Rule (h = ',num2str(silvermansBW,'%.3f'),')'],...
    ['Critical Bandwidth (h* = ',num2str(criticalBW,'%.3f'),')']},...
    'box','off','location','northwest')
xlim(prctile(grid,[0 100])); ylim([0 0.07])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 3 1.25])
xlabel('projection onto PC1','FontSize',8)
ylabel('PDF','FontSize',8)
exportgraphics(gcf,[saveDir,'figS11b.pdf'],'ContentType','vector')

%% Figure S11C - Significance of Silverman's critical bandwidth

% Bootstrap resample data to get significance of critical bandwidth
rng(1);
numReps = 2000;
criticalBoot = nan(numReps,1);
numPoints = length(PC1);
for i = 1:numReps
    % Resample
    yResample = datasample(PC1,numPoints);

    % Rescale
    noise = randn(numPoints,1);
    yRescale = (yResample + criticalBW*noise/sqrt(1 + criticalBW^2/var(PC1)));

    % Find new critical value for each resample
    [criticalBoot(i)] = silvermansTest(yRescale,false);
end

% Calculate p-value
pValue = sum([criticalBoot;criticalBW] >= criticalBW)/(numReps+1)

figure; hold on
histogram(criticalBoot,36,...
    'Normalization','probability','FaceColor',colorValue.*[1 1 1])
xline(criticalBW,'LineWidth',2,'Color',exploreColor,'Alpha',1)
ylim([0 0.1]); xlim([0.1 0.42])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 3 1.25])
xlabel('critical bandwidth (h*) of bootstrap resamples','FontSize',8)
ylabel('probability','FontSize',8)
text(criticalBW-0.01,0.085,{['h* = ',num2str(criticalBW,'%.3f')];['p < ',num2str(pValue,'%.3f')]},...
    'FontSize',8,'Color',exploreColor,'HorizontalAlignment','right','FontName','Arial')
exportgraphics(gcf,[saveDir,'figS11c.pdf'],'ContentType','vector')

%% Figure S11D - Explore vs. Exploit Clusters on all data

exploitPosterior = encounter.exploitPosterior(include);
sensePosterior = encounter.sensePosterior(include);

% Get contours of standard deviation(s) of gaussians
alphaVal = [0.4 0.3 0.2];
stdValues = [1 2 3];
numPoints = 1000;
exploitEllipses = gaussianContours(exploitDist.mu,exploitDist.Sigma,stdValues,numPoints);
exploreEllipses = gaussianContours(exploreDist.mu,exploreDist.Sigma,stdValues,numPoints);

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
for i = 1:length(alphaVal)
    patch(exploitEllipses(:,2*i-1),exploitEllipses(:,2*i),...
        min(exploitColor + (1 - exploitColor).*(colorValue-0.1),1),...
        'FaceAlpha',alphaVal(i),'EdgeAlpha',0)
    patch(exploreEllipses(:,2*i-1),exploreEllipses(:,2*i),...
        min(sampleColor + (1 - sampleColor).*(colorValue-0.1),1),...
        'FaceAlpha',alphaVal(i),'EdgeAlpha',0)
end
scatter(clusterVars(:,1),clusterVars(:,2),1,'CData',...
    exploitPosterior.*exploitColor + (1-exploitPosterior).*sampleColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none')
xlim([-2.1 2.1]); ylim([1.2 2.8]); legend('off')
xticks(log10([2 10 60 600 3600]./60)); xticklabels({'1/30','1/6','1','10','60'});
yticks(log10([25 50 100 200 400])); yticklabels([25 50 100 200 400])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[1 1 3 3])
xlabel('encounter duration (min)','FontSize',8); ylabel('avg. velocity on patch (μm/s)','FontSize',8)

colors = interp1([0;1],[sampleColor;exploitColor],linspace(0,1,100));
colormap(colors); c = colorbar('Units','inches');
set(c,'Position',[4.1 1 0.15 3]);
set(c.Label,'String','p(exploit cluster)','FontSize',8,'Rotation',-90)
exportgraphics(gcf,[saveDir,'figS11d.pdf'],'ContentType','vector')

% [X,Y] = meshgrid(linspace(-2.1,2.1,1000),linspace(1.2,2.8,1000));
% posteriorGrid = posterior(exploitModel,[X(:),Y(:)]);
% exploitPDF = reshape(posteriorGrid(:,exploitCluster),size(X));
% figure; imagesc(exploitPDF); colormap(colors)

%% Figure S11E - Posterior variance vs. alpha

posteriorVar = reshape(permute(posteriorVarExploit,[2,1,3]),length(alpha),[]);
posteriorVarCI = prctile(posteriorVar,[2.5 50 97.5],2);

figure; hold on;
patch([alpha,flip(alpha)],[posteriorVarCI(:,1)',flip(posteriorVarCI(:,3))'],...
    [0.5 0.5 0.5],'EdgeColor','none')
plot(alpha,posteriorVarCI(:,2),'k','LineWidth',1)
scatter(alpha(bestModel),posteriorVarCI(bestModel,2),'k','filled')
text(alpha(bestModel),posteriorVarCI(bestModel,2)+0.003,['α = ',num2str(alpha(bestModel))],...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',8,'FontName','Arial')
xlim(prctile(alpha,[0 100]))
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
ax = gca; set(gca,'Position',[ax.Position(1:2) 1 2])
xlabel('α','FontSize',8); ylabel('posterior variance','FontSize',8)
exportgraphics(gca,[saveDir,'figS11e.pdf'])

%% Figure S11F - Sorted posterior probability of p(exploit) for all data

figure; hold on;
plot(sort(1-exploitPosterior,'descend'),'LineWidth',2,'Color',sampleColor)
plot(sort(exploitPosterior,'ascend'),'LineWidth',2,'Color',exploitColor)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlim([0 length(exploitPosterior)]); xticks([0 length(exploitPosterior)]);
xticklabels({'0',num2str(length(exploitPosterior))}); yticks([0 1]);
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('# encounters','FontSize',8); ylabel('p(cluster)','FontSize',8)
exportgraphics(gca,[saveDir,'figS11f.pdf'])

%% Figure S11G - Sorted posterior probability of p(exploit) for data in Fig 2

exploitPosteriorConc = encounter.exploitPosterior(include & indEncounter);

figure; hold on;
plot(sort(1-exploitPosteriorConc,'descend'),'LineWidth',2,'Color',sampleColor)
plot(sort(exploitPosteriorConc,'ascend'),'LineWidth',2,'Color',exploitColor)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlim([0 length(exploitPosteriorConc)]); xticks([0 length(exploitPosteriorConc)]);
xticklabels({'0',num2str(length(exploitPosteriorConc))}); yticks([0 1]);
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('# encounters','FontSize',8); ylabel('p(cluster)','FontSize',8)
exportgraphics(gca,[saveDir,'figS11g.pdf'])

%% Figure S12-1A - Quantification of Deceleration, Min. Velocity, Max Change

wormNum = 360;
encounterNum = 6;
exampleEncounter = find(indEncounter & encounter.wormNum == 360 & encounter.id == encounterNum);
indData = data.wormNum == wormNum;
time = data.time(indData);
velocitySmooth = data.velocitySmooth(indData);

% Compute deceleration at slowdown via linear fit from -1.5 to +6.5
% seconds from patch entry
indDecelerate = time >= encounter.timeEnter(exampleEncounter) - 1.5 & ...
    time <= encounter.timeEnter(exampleEncounter) + 6.5;
dT = time(indDecelerate) - mean(time(indDecelerate),'omitnan');
dV = velocitySmooth(indDecelerate) - mean(velocitySmooth(indDecelerate),'omitnan');
decelerate = sum(dT.*dV,'omitnan')./sum(dT.^2,'omitnan');

% Smooth the velocity using cubic spline interpolation and compute gradient
velocitySpline = csaps(time,velocitySmooth,0.5,time);
accelerationSpline = gradient(velocitySpline,time);

% Find peak before
peakBefore = find(accelerationSpline >= 0 & time <= encounter.timeEnter(exampleEncounter) & ...
    time >= encounter.timeEnter(exampleEncounter) - 10,1,'last');
peakBefore = find(abs(time - time(peakBefore)) <= 1);
[velocityBeforeEnter,indBefore] = max(velocitySmooth(peakBefore));

% Find min velocity on
indOn = find(time >= encounter.timeEnter(exampleEncounter) & ...
    time <= encounter.timeExit(exampleEncounter));
onPatchVelocity = velocitySmooth(indOn);
[minVelocity,indMin] = min(onPatchVelocity);

figure; hold on;

% Plot encounter & velocity
patch(reshape(repmat([encounter.timeEnter(exampleEncounter) encounter.timeExit(exampleEncounter)],2,1),4,1),...
    [0 400 400 0]',colorValue.*[1 1 1],'EdgeAlpha',0,'FaceAlpha',0.5)
plot(time,velocitySmooth,'k','LineWidth',1)
xlim([encounter.timeEnter(exampleEncounter)-20 encounter.timeExit(exampleEncounter)+10])
xticks(encounter.timeEnter(exampleEncounter) + [-20:20:100]);
xticklabels(-20:20:100); ylim([0 400])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
xlabel('time from patch entry (s)','FontSize',8); ylabel('velocity (μm/s)','FontSize',8)
ax = gca; set(gca,'Position',[ax.Position(1:2) 3.5 1.5])

% Show deceleration
xline(time(find(indDecelerate,1,'first')),'Color','#77AC30')
xline(time(find(indDecelerate,1,'last')),'Color','#A2142F')
f = fit(time(indDecelerate),velocitySmooth(indDecelerate),'poly1');
plot(time(indDecelerate),f(time(indDecelerate)),'Color','#7E2F8E','LineWidth',1)
text(mean(time(indDecelerate),'omitnan'),mean(velocitySmooth(indDecelerate)+30,'omitnan'),...
    'Δ velocity (μm/s^2)','HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Color','#7E2F8E','FontName','Arial','FontSize',8);
text(time(find(indDecelerate,1,'first'))-1,360,...
    '-1.5','HorizontalAlignment','right','VerticalAlignment','bottom',...
    'Color','#77AC30','FontName','Arial','FontSize',8);
text(time(find(indDecelerate,1,'last'))+1,360,...
    '+6.5','HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Color','#A2142F','FontName','Arial','FontSize',8);

% Show min velocity & max before
scatter(time(peakBefore(indBefore)),velocityBeforeEnter,'MarkerEdgeColor','#A2142F')
text(time(peakBefore(indBefore))+2,velocityBeforeEnter,...
    'max. velocity before encounter (μm/s)','HorizontalAlignment','left','VerticalAlignment','middle',...
    'Color','#A2142F','FontName','Arial','FontSize',8);
scatter(time(indOn(indMin)),minVelocity,'s','MarkerEdgeColor','#0072BD')
text(time(indOn(indMin))+2,minVelocity-10,...
    'min. velocity on patch (μm/s)','HorizontalAlignment','left','VerticalAlignment','middle',...
    'Color','#0072BD','FontName','Arial','FontSize',8);
plot([200 200],[minVelocity velocityBeforeEnter],'Color','#D95319','LineWidth',1)
plot([199 201],[minVelocity minVelocity],'Color','#D95319','LineWidth',1)
plot([199 201],[velocityBeforeEnter velocityBeforeEnter],'Color','#D95319','LineWidth',1)
text(196,mean([minVelocity velocityBeforeEnter]),sprintf('max. Δ velocity (μm/s)'),...
    'Color','#D95319','FontSize',8,'FontName','Arial','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90)

exportgraphics(gcf,[saveDir,'figS12-1a.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S12-1B - Deceleration

figure; hold on
xValues = 1:numGroups; xValues(end-1:end) = xValues(end-1:end) + 0.5;
yline(0,'k')
for i = 1:numGroups
    ind = indEncounter & G == conditionG(i) & ~minVelocityOnly;
    deceleration = encounter.decelerate(ind); %encounter.velocityBeforeEnter(ind)-encounter.velocityOnMin(ind);
    v = Violin({deceleration},xValues(i),'ViolinColor',{colorValues(i).*[1 1 1]},...
        'QuartileStyle','shadow','ShowMean',true,'HalfViolin','full');
    v.ScatterPlot.SizeData = 5; v.MedianPlot.SizeData = 10;
end
xlim([0.25 numGroups+1.25])
xticks(xValues); xticklabels({'0','0.05','0.1','0.5','1','2',...
    '3','4','5','10','1','1'})
%yticks(0:150:450); ylim([0 450])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 2 1])
ylabel('Δ velocity (μm/s^2)','FontSize',8)
exportgraphics(gcf,[saveDir,'figS12-1b.pdf'],'BackgroundColor','none','ContentType','vector')

%% Figure S12-1C - Heatmap of Min. Velocity vs. Deceleration

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
h1 = histogram2(encounter.velocityOnMin(indEncounter),encounter.decelerate(indEncounter),...
    linspace(0,450,51),linspace(-65,45,51),'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','probability');
h2 = histcounts2(encounter.velocityOnMin(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    linspace(0,450,51),linspace(0,450,51),'Normalization','probability');
legend('off'); xlim([0 450]); ylim([-65 45]); xticks(0:100:400); yticks(-60:20:40)
maxVal = max([h1.Values,h2],[],'all');
clim([0 maxVal]); colormap(flipud(gray))
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[1 1 3 3])
xlabel('min. velocity on patch (μm/s)','FontSize',8,'Color','#0072BD');
ylabel('Δ velocity (μm/s^2)','FontSize',8,'Color','#7E2F8E')
exportgraphics(gcf,[saveDir,'figS12-1c.pdf'],'ContentType','vector')

%% Figure S12-1D - Heatmap of Min. Velocity vs. Max Change

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
histogram2(encounter.velocityOnMin(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    linspace(0,450,51),linspace(0,450,51),'DisplayStyle','tile','ShowEmptyBins','on',...
    'Normalization','probability');
legend('off'); xlim([0 450]); ylim([0 450]); xticks(0:100:400); yticks(0:100:400)
clim([0 maxVal]); colormap(flipud(gray))
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[1 1 3 3])
xlabel('min. velocity on patch (μm/s)','FontSize',8,'Color','#0072BD');
ylabel('max. Δ velocity (μm/s)','FontSize',8,'Color','#D95319')

c = colorbar('Units','inches');
set(c,'Position',[4.1 1 0.15 3]);
set(c.Label,'String','p(encounter)','FontSize',8,'Rotation',-90)
exportgraphics(gcf,[saveDir,'figS12-1d.pdf'],'ContentType','vector')

%% Figure S12-1E - Projection onto PC1 and Silverman's Critical Bandwidth

% Create mean-centered matrix of variables A
A = [encounter.velocityOnMin(indEncounter),encounter.decelerate(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter)];
A = A(~any(isnan(A),2),:);
A = A - mean(A,1);

% Compute first eigenvectors of A
[V,~] = eigs(A'*A,1);

% Project values in A onto the eigenvector V
PC1 = A*V(:);

% Get critical bandwidth of data
[criticalBW,silvermansBW,grid] = silvermansTest(PC1,false);

% Plot histogram of values projected onto PC1
figure; hold on;
histogram(PC1,'BinEdges',grid,'Normalization','probability','FaceColor',colorValue.*[1 1 1])
binWidth = grid(2) - grid(1);
p1 = plot(grid,binWidth*mvksdensity(PC1,grid,'Bandwidth',silvermansBW),'LineWidth',2,'Color',exploitColor);
p2 = plot(grid,binWidth*mvksdensity(PC1,grid,'Bandwidth',criticalBW),'LineWidth',2,'Color',exploreColor);
legend([p1,p2],{['Silverman''s Rule (h = ',num2str(silvermansBW,'%.3f'),')'],...
    ['Critical Bandwidth (h* = ',num2str(criticalBW,'%.3f'),')']},...
    'box','off','location','northeast')
xlim(prctile(grid,[0 100]));
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 3 1.25])
xlabel('projection onto PC1','FontSize',8)
ylabel('PDF','FontSize',8)
exportgraphics(gcf,[saveDir,'figS12-1e.pdf'],'ContentType','vector')

%%
% [V,D] = eigs(A'*A,2);
% 
% W = null(V);
% 
% [X2,Y2] = meshgrid(-200:10:400,-50:5:50);
% Z2=(1/V(3,1)).*((V(1,1).*X2)+(V(2,1).*Y2));
% 
% figure; hold on
% scatter3(A(:,1),A(:,2),A(:,3),'k')
% %surf(X2,Y2,Z2)
% D2 = sqrt(diag(D));
% plot3([0,V(1,1)*D2(1)],[0,V(2,1)*D2(1)],[0,V(3,1)*D2(1)],'r','LineWidth',3)
% plot3([0,V(1,2)*D2(2)],[0,V(2,2)*D2(2)],[0,V(3,2)*D2(2)],'b','LineWidth',3)
% 
% k = sum(-V(:,1)'.*A,2)./sum(V(:,1).^2);
% A2 = A + k*V(:,1)';
% 
% figure; hold on;
% scatter3(A(:,1),A(:,2),A(:,3),'k')
% scatter3(A2(:,1),A2(:,2),A2(:,3),'b')
% plot3([A(:,1),A2(:,1)]',[A(:,2),A2(:,2)]',[A(:,3),A2(:,3)]','k','LineWidth',0.1)
% surf(X2,Y2,Z2)

%% Figure S12-1F - Significance of Silverman's critical bandwidth

% Bootstrap resample data to get significance of critical bandwidth
rng(1);
numReps = 2000;
criticalBoot = nan(numReps,1);
numPoints = length(PC1);
for i = 1:numReps
    % Resample
    yResample = datasample(PC1,numPoints);

    % Rescale
    noise = randn(numPoints,1);
    yRescale = (yResample + criticalBW*noise/sqrt(1 + criticalBW^2/var(PC1)));

    % Find new critical value for each resample
    [criticalBoot(i)] = silvermansTest(yRescale,false);
end

% Calculate p-value
pValue = sum([criticalBoot;criticalBW] >= criticalBW)/(numReps+1)

figure; hold on
histogram(criticalBoot,'BinEdges',10:1:50,...
    'Normalization','probability','FaceColor',colorValue.*[1 1 1])
xline(criticalBW,'LineWidth',2,'Color',exploreColor,'Alpha',1)
ylim([0 0.07]); xlim([10 50])
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[ax.Position(1:2)+1 3 1.25])
xlabel('critical bandwidth (h*) of bootstrap resamples','FontSize',8)
ylabel('probability','FontSize',8)
text(criticalBW-1,0.06,{['h* = ',num2str(criticalBW,'%.3f')];['p = ',num2str(pValue,'%.3f')]},...
    'FontSize',8,'Color',exploreColor,'HorizontalAlignment','right','FontName','Arial')
exportgraphics(gcf,[saveDir,'figS12-1f.pdf'],'ContentType','vector')

%% Figure S12-2A - Sample vs. Search Clusters for all data

% Get QDA paraboloid boundary
K = senseModel.Coeffs(1,2).Const;
L = senseModel.Coeffs(1,2).Linear; 
Q = senseModel.Coeffs(1,2).Quadratic;
quadratic3D = @(x1,x2,x3) K + L(1)*x1 + L(2)*x2 + L(3)*x3 + ...
    Q(1,1)*x1.^2 + Q(2,2)*x2.^2 + Q(3,3)*x3.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + (Q(1,3)+Q(3,1))*x1.*x3 + (Q(2,3)+Q(3,2))*x2.*x3;
x2_median = prctile(encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),50);
x3_median = prctile(encounter.decelerate(indEncounter),50);
quadratic2D_13_median = @(x1,x3) K + L(1)*x1 + L(2)*x2_median + L(3)*x3 + ...
    Q(1,1)*x1.^2 + Q(2,2)*x2_median.^2 + Q(3,3)*x3.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2_median + (Q(1,3)+Q(3,1))*x1.*x3 + (Q(2,3)+Q(3,2))*x2_median.*x3;
quadratic2D_12_median = @(x1,x2) K + L(1)*x1 + L(2)*x2 + L(3)*x3_median + ...
    Q(1,1)*x1.^2 + Q(2,2)*x2.^2 + Q(3,3)*x3_median.^2 + ...
    (Q(1,2)+Q(2,1))*x1.*x2 + (Q(1,3)+Q(3,1))*x1.*x3_median + (Q(2,3)+Q(3,2))*x2.*x3_median;

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
fp = fimplicit(quadratic2D_13_median,[0 450 -65 45],'LineStyle','none');
patch([450 fp.XData 450],[-65 fp.YData 45],min(searchColor + (1 - searchColor).*(colorValue-0.1),1),...
         'FaceAlpha',min(alphaVal),'EdgeAlpha',0)
patch([0 fp.XData 0],[-65 fp.YData 45],min(sampleColor + (1 - sampleColor).*(colorValue-0.1),1),...
         'FaceAlpha',min(alphaVal),'EdgeAlpha',0)
scatter(encounter.velocityOnMin(indEncounter),encounter.decelerate(indEncounter),...
    2,'CData',encounter.sensePosterior(indEncounter).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none')
scatter(encounter.velocityOnMin(indEncounter & indSense),encounter.decelerate(indEncounter & indSense),...
    2,'CData',encounter.sensePosterior(indEncounter & indSense).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter & indSense)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor',exploitColor,'LineWidth',0.5)
scatter(encounter.velocityOnMin(indEncounter & indNonsense),encounter.decelerate(indEncounter & indNonsense),...
    2,'CData',encounter.sensePosterior(indEncounter & indNonsense).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter & indNonsense)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor',exploreColor,'LineWidth',0.5)
legend('off'); xlim([0 450]); ylim([-65 45]); xticks(0:100:400); yticks(-60:20:40)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[1 1 3 3])
xlabel('min. velocity on patch (μm/s)','FontSize',8,'Color','#0072BD');
ylabel('Δ velocity (μm/s^2)','FontSize',8,'Color','#7E2F8E')
exportgraphics(gcf,[saveDir,'figS12-2a.pdf'],'ContentType','vector')

%% Figure S12-2B - Sample vs. Search Clusters for all data (alternate view)

figure('Position',[400 200 560 560]);
subplot(6,6,[1:3,7:9,13:15]); hold on
fp = fimplicit(quadratic2D_12_median,[0 450 0 450],'LineStyle','none');
patch([450 fp.XData 450],[0 fp.YData 450],min(searchColor + (1 - searchColor).*(colorValue-0.1),1),...
         'FaceAlpha',min(alphaVal),'EdgeAlpha',0)
patch([0 fp.XData 0],[0 fp.YData 450],min(sampleColor + (1 - sampleColor).*(colorValue-0.1),1),...
         'FaceAlpha',min(alphaVal),'EdgeAlpha',0)
scatter(encounter.velocityOnMin(indEncounter),encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    1,'CData',encounter.sensePosterior(indEncounter).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none')
scatter(encounter.velocityOnMin(indEncounter & indSense),...
    encounter.velocityBeforeEnter(indEncounter & indSense)-encounter.velocityOnMin(indEncounter & indSense),...
    2,'CData',encounter.sensePosterior(indEncounter & indSense).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter & indSense)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor',exploitColor,'LineWidth',0.5)
scatter(encounter.velocityOnMin(indEncounter & indNonsense),...
    encounter.velocityBeforeEnter(indEncounter & indNonsense)-encounter.velocityOnMin(indEncounter & indNonsense),...
    2,'CData',encounter.sensePosterior(indEncounter & indNonsense).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter & indNonsense)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor',exploreColor,'LineWidth',0.5)
legend('off'); xlim([0 450]); ylim([0 450]); xticks(0:100:400); yticks(0:100:450)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out');
ax = gca; set(gca,'Position',[1 1 3 3])
xlabel('min. velocity on patch (μm/s)','FontSize',8,'Color','#0072BD');
ylabel('max. Δ velocity (μm/s)','FontSize',8,'Color','#D95319')

colors = interp1([0;1],[searchColor;sampleColor],linspace(0,1,100));
colormap(colors); c = colorbar('Units','inches');
set(c,'Position',[4.1 1 0.15 3]);
set(c.Label,'String','p(sense cluster)','FontSize',8,'Rotation',-90)
exportgraphics(gcf,[saveDir,'figS12-2b.pdf'],'ContentType','vector')

%% Figure S12-2C - 3D Sense Clustering

figure('Position',[400 200 560 560]); hold on
fp = fimplicit3(quadratic3D,[0 450 0 450 -65 45],'LineStyle','none','FaceColor','k','FaceAlpha',0.2);
scatter3(encounter.velocityOnMin(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    encounter.decelerate(indEncounter),5,...
    'CData',encounter.sensePosterior(indEncounter).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
legend('off'); xlim([0 450]); ylim([0 450]); zlim([-65 45]);
xticks(0:100:400); yticks(0:100:450); zticks(-60:20:40)
xlabel('min. velocity on patch (μm/s)','FontSize',8,'Color','#0072BD');
ylabel('max. Δ velocity (μm/s)','FontSize',8,'Color','#D95319')
zlabel('Δ velocity (μm/s^2)','FontSize',8,'Color','#7E2F8E')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out',...
    'XGrid','on','YGrid','on','ZGrid','on','Color','w','Projection','perspective');
ax = gca; set(gca,'Position',[1 1 3 3])
view([-30 10])
exportgraphics(gcf,[saveDir,'figS12-2c.pdf'],'ContentType','vector')

%% Video S5 - 3D animation of QDA for Sense-Nonsense

figure('Position',[100 100 1080 1080],'Color','k'); hold on
fp = fimplicit3(quadratic3D,[0 450 0 450 -65 45],'LineStyle','none','FaceColor','w','FaceAlpha',0.2);
scatter3(encounter.velocityOnMin(indEncounter),...
    encounter.velocityBeforeEnter(indEncounter)-encounter.velocityOnMin(indEncounter),...
    encounter.decelerate(indEncounter),25,...
    'CData',encounter.sensePosterior(indEncounter).*sampleColor + ...
    (1-encounter.sensePosterior(indEncounter)).*searchColor,...
    'MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% scatter3(encounter.velocityOnMin(indEncounter & indSense),...
%     encounter.velocityBeforeEnter(indEncounter & indSense)-encounter.velocityOnMin(indEncounter & indSense),...
%     2,'CData',exploitColor,...
%     'MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
% scatter3(encounter.velocityOnMin(indEncounter & indNonsense),...
%     encounter.velocityBeforeEnter(indEncounter & indNonsense)-encounter.velocityOnMin(indEncounter & indNonsense),...
%     2,'CData',exploreColor,...
%     'MarkerFaceColor','flat','MarkerEdgeColor','none','MarkerFaceAlpha',0.2)
legend('off'); xlim([0 450]); ylim([0 450]); zlim([-65 45]);
xticks(0:100:400); yticks(0:100:450); zticks(-60:20:40)
xlabel('min. velocity on patch (μm/s)');
ylabel('max. Δ velocity (μm/s)')
zlabel('Δ velocity (μm/s^2)')
set(gca,'FontName','Arial','FontSize',18,'TickDir','out',...
    'Color','k','XColor','w','YColor','w','ZColor','w','BoxStyle','full');
set(gcf,'ToolBar','none')

videoName = fullfile(path,'figures','VideoSupp','VideoS5');
v = VideoWriter(videoName,'MPEG-4');
v.FrameRate = 15; v.Quality = 100;
open(v);

move1 = interp1([0;1],[0 0;-30 10],linspace(0,1,30));
move2 = interp1([0;1],[-30 10;-390 10],linspace(0,1,195));
move3 = interp1([0;1],[-30 10;0 90],linspace(0,1,45));
%move4 = interp1([0;1],[0 90;0 0],linspace(0,1,45));
pause1 = repmat([-30 10],15,1);
%pause2 =  repmat([0 90],15,1);
moves = [move1;pause1;move2;pause1;move3];
for i = 1:length(moves)
    view(moves(i,:))
    frame = getframe(gcf);
    writeVideo(v,frame);
end

% Close video
writeVideo(v,frame);
close(v)

%% Figure S12-2D - Sorted posterior probability of p(sense) for all data

figure; hold on;
plot(sort(1-sensePosterior,'descend'),'LineWidth',2,'Color',searchColor)
plot(sort(sensePosterior,'ascend'),'LineWidth',2,'Color',sampleColor)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlim([0 length(sensePosterior)]); xticks([0 length(sensePosterior)]);
xticklabels({'0',num2str(length(sensePosterior))}); yticks([0 1]);
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('# encounters','FontSize',8); ylabel('p(cluster)','FontSize',8)
exportgraphics(gca,[saveDir,'figS12-2d.pdf'])

%% Figure S12-2E - Sorted posterior probability of p(sense) for data in Fig 2

sensePosteriorConc = encounter.sensePosterior(include & indEncounter);

figure; hold on;
plot(sort(1-sensePosteriorConc,'descend'),'LineWidth',2,'Color',searchColor)
plot(sort(sensePosteriorConc,'ascend'),'LineWidth',2,'Color',sampleColor)
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out'); 
xlim([0 length(sensePosteriorConc)]); xticks([0 length(sensePosteriorConc)]);
xticklabels({'0',num2str(length(sensePosteriorConc))}); yticks([0 1]);
ax = gca; set(gca,'Position',[ax.Position(1:2) 1.5 0.75])
xlabel('# encounters','FontSize',8); ylabel('p(cluster)','FontSize',8)
exportgraphics(gca,[saveDir,'figS12-2e.pdf'])

%% Figure S12-2F -  Estimate of posterior probabilities from min. velocity only via marginalization

% figure; hold on
% scatter(clusterVar1,posteriorProbSense(minVelocityOnly,senseCluster),5,...
%     'CData',posteriorProbSenseMarginal.*exploreColor + ...
%     (1-posteriorProbSenseMarginal).*nonSenseColor,...
%     'MarkerFaceColor','flat','MarkerEdgeColor','none')
% xlim([0 350])
% set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out')
% xlabel('min. velocity on patch (μm/s)','FontSize',8); ylabel('censored p(sense)','FontSize',8)
% ax = gca; set(gca,'Position',[1 1 1.5 1.5])
% 
% colors = interp1([0;1],[nonSenseColor;exploreColor],linspace(0,1,100));
% colormap(colors); c = colorbar('Units','inches');
% set(c,'Position',[2.6 1 0.15 1.5]);
% set(c.Label,'String','marginalized p(sense)','FontSize',8,'Rotation',-90)
% exportgraphics(gca,[saveDir,'figS11j.pdf'])

figure; hold on
patch([0 1 1 0],[0.5 0.5 1 1],sampleColor,'FaceAlpha',min(alphaVal),'EdgeColor','none')
patch([0 1 1 0],[0.5 0.5 0 0],searchColor,'FaceAlpha',min(alphaVal),'EdgeColor','none')
xline(0.5,'k--'); yline(0.5,'k--')
scatter(posteriorProbSense(minVelocityOnly,senseCluster),posteriorProbSenseMarginal,...
    20,'CData',clusterVar1,'MarkerFaceColor','flat','MarkerEdgeColor','none')
set(gca,'Units','inches','FontName','Arial','FontSize',8,'Box','on','TickDir','out',...
    'XTick',0:0.5:1,'YTick',0:0.5:1)
xlabel('censored p(sense)','FontSize',8); ylabel('marginalized p(sense)','FontSize',8)
ax = gca; set(gca,'Position',[1 1 1.5 1.5]);

colormap(flipud(jet)); c = colorbar('Units','inches');
set(c,'Position',[2.6 1 0.15 1.5],'Limits',[0 300],'Direction','reverse',...
    'Ticks',0:100:300,'TickLabels',300:-100:0);
set(c.Label,'String','min. velocity on patch (μm/s)','FontSize',8,'Rotation',-90)
exportgraphics(gca,[saveDir,'figS12-2f.pdf'])