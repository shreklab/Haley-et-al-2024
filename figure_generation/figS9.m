% Load data

path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
expName = 'foragingConcentration';
bodyPart = 'midpoint';
load(fullfile(path,expName,'experimentInfo.mat'),'info');
load(fullfile(path,expName,[bodyPart,'.mat']),'data');
saveDir = [path,'figures\FigS9\'];

%% Figure S9 - Example traces

wormNum = [465, 474, 467, 718, 491, 695, 478, 343, 487, 360, 662, 624]; % [0, 0.05, 0.1, 0.5, 1, 2, 4, 5, 10, 1(48H)]
for i = 1:length(wormNum)
    % Edit code to save pngs & increase lineWidth to 3
    plotTracks(info,data,wormNum(i),'timeOffset',saveDir,bodyPart,'yes');
end