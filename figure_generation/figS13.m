%% Load Data

expName = 'foragingConcentration';
path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');

saveDir = [path,'figures\FigS13\'];

%% Figure S13 - Example traces (Large Single Patches)
wormNum = [515,546,529,533,514,509,526,519,539,537,506]; % 0,0.05,0.1,0.5,1,2,3,4,5,10,1(48H)
for i = 1:length(wormNum)
    % Edit code to save pngs & increase lineWidth to 3
    plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
end