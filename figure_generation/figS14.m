%% Load Data

expName = 'foragingMini';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');

saveDir = [path,'figures\FigS14\'];

%% Figure S14 - Example traces (Large Single Patches)
wormNum = [212,211,216,221,214,218]; % 0,0.5,1,5,10,1(48H)
for i = 1:length(wormNum)
    % Edit code to save pngs & increase lineWidth to 3
    plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
end