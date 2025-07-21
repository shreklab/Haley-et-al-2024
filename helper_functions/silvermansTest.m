function [criticalBW,silvermansBW,grid] = silvermansTest(data,plotData,gridSize,tolerance)
% [CRITICALBW,SILVERMANSBW,GRID] = SILVERMANSTEST(DATA,PLOTDATA,GRIDSIZE,TOLERANCE)
%
%   SILVERMANSTEST calculates two estimates for the optimal kernel density
%   estimation (KDE) bandwidth for a given dataset. It computes a
%   bandwidth using Silverman's rule of thumb and iteratively finds the
%   "critical bandwidth" where the number of modes (peaks) in the KDE
%   switches from multimodal to unimodal.
%
%   INPUTS:
%       - data [NxD double]: The input data matrix, with N observations
%           and D dimensions.
%       - plotData [logical]: A flag to control whether a plot is
%           generated comparing the results (true = plot, false = no plot).
%       - gridSize [double]: (Optional) The number of points in the grid
%           used for KDE evaluation. Defaults to 101.
%       - tolerance [double]: (Optional) The desired precision for the
%           iterative critical bandwidth calculation. Defaults to 1e-3.
%
%   OUTPUTS:
%       - criticalBW [double]: The calculated critical bandwidth where the
%           KDE becomes unimodal.
%       - silvermansBW [double]: The bandwidth calculated using Silverman's
%           rule of thumb.
%       - grid [1xG double]: The grid of points over which the KDE was
%           evaluated, where G is the gridSize.
%
%   Written 7/11/2025 by Jess Haley in MATLAB R2024a.
%
%   See also MVKSDENSITY, FINDPEAKS, PRCTILE, HISTOGRAM.

% Get range of data and round outwards
minMax = prctile(data,[0 100]);
orderMag = floor(log10(abs(minMax)));
dataRange = minMax.*(10.^-orderMag);
dataRange(1) = floor(dataRange(1)); dataRange(2) = ceil(dataRange(2));
dataRange = dataRange.*(10.^orderMag);

% Estimate a good bandwidth using Silverman's rule
d = size(data,2); % # dimensions
n = size(data,1); % # observations
silvermansBW = std(data,[],1,'omitnan').*(4/(n*(d+2))).^(1/(d+4)); % Silverman's rule

% Calculate critical bandwidth where the KDE mode switches from 1 to 2
if nargin < 3
    gridSize = 101;
end
if nargin < 4
    tolerance = 1e-3;
end
precision = tolerance*10;
numBW = 20;
rangeBW = silvermansBW*[1/10,10]; % initialize range of bandwidths to try
grid = linspace(dataRange(1),dataRange(2),gridSize); % range of data

while precision > tolerance
    bw = linspace(rangeBW(1),rangeBW(2),numBW);
    numPeaks = nan(numBW,1);

    % Calculate # of peaks (modes) for each bw
    for i = 1:numBW
        f = mvksdensity(data,grid,'Bandwidth',bw(i));
        [pks,~] = findpeaks(f);
        numPeaks(i) = numel(pks);
        if numPeaks(i) == 1
            break
        end
    end

    % Find critical bandwidth where the # of modes goes down to 1
    criticalInd = find(-diff(numPeaks==1).*diff(numPeaks>1));
    rangeBW = [bw(criticalInd) bw(criticalInd+1)];

    % Calculate precision of critical bandwidth
    precision = rangeBW(2)-rangeBW(1);
end
criticalBW = mean(rangeBW);

% Plot histogram and KDEs of data
if plotData
    figure; hold on;
    histogram(data,'BinEdges',grid,'Normalization','probability')
    binWidth = grid(2) - grid(1);
    p1 = plot(grid,binWidth*mvksdensity(data,grid,'Bandwidth',silvermansBW),'LineWidth',2);
    p2 = plot(grid,binWidth*mvksdensity(data,grid,'Bandwidth',criticalBW),'LineWidth',2);
    legend([p1,p2],{['Silvermans Rule (h = ',num2str(silvermansBW),')'],...
        ['Critical Bandwidth (h = ',num2str(criticalBW),')']})
    xlabel('Projection onto PC1')
    ylabel('PDF')
end

end