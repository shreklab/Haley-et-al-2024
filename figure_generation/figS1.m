%% Load Data

path = 'Z:\jhaley\foragingPaper\';
addpath(genpath(path))
saveDir = [path,'figures\FigS1\'];
load([path,'foragingConcentration\videos\22-03-18\data\lawn\2022-03-18_12-55-44_2.mat']);

%% Figure S1A - Biorender Graphic

% Created graphic in Biorender

%% Figure S1B - Contrast Frame(s)

imwrite(read(lawn.video.vidObject,100),[saveDir,'figS1b_1.png'])
imwrite(read(lawn.video.vidObject,101),[saveDir,'figS1b_2.png'])
imwrite(read(lawn.video.vidObject,102),[saveDir,'figS1b_3.png'])
imwrite(read(lawn.video.vidObject,103),[saveDir,'figS1b_4.png'])
imwrite(read(lawn.video.vidObject,104),[saveDir,'figS1b_5.png'])
imwrite(read(lawn.video.vidObject,105),[saveDir,'figS1b_6.png'])

%% Figure S1C - Filter Image

thresh = 40;
filterImage = min(lawn.filter.image,thresh);
filterImage = imcomplement(uint8(filterImage*256/thresh));
imwrite(filterImage,[saveDir,'figS1c.png'])

%% Figure S1D - Patch & Arena Location

% Used same image as in Fig 1A - slight offset between Fig S1C due to image
% registration that was used

%% Fig S1E - Behavior Frame(s)

% Downloaded frame 7000 of video of worm in Fig 1A using WormLab

%% Fig S1F - Behavior Tracking

% Downloaded annotated frame 7000 of video of worm in Fig 1A using WormLab

%% Fig S1G - Midbody location over time

% Same image as in Fig 1A