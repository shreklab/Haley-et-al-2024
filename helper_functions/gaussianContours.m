function [ellipses] = gaussianContours(mu,sigma,alpha,numPoints)
% [ellipses] = GAUSSIANCONTOURS(mu,sigma,alpha,numPoints)
%
%   GAUSSIANCONTOURS calculates the coordinates for plotting concentric
%   ellipses representing the standard deviation contours of a 2D Gaussian
%   distribution.
%
%   INPUTS:
%       - mu [1x2 double]: A vector containing the mean [X,Y] of the
%           Gaussian distribution.
%       - sigma [2x2 double]: The covariance matrix of the Gaussian
%           distribution.
%       - alpha [1xN double]: An array of standard deviation multipliers
%           to define the contour levels (e.g., [1, 2, 3] for 1, 2, and 3
%           standard deviations).
%       - numPoints [double]: The number of points to use for drawing
%           each ellipse contour.
%
%   OUTPUTS:
%       - ellipses [numPoints x 2*N double]: A matrix where each pair of
%           columns [X, Y] contains the coordinates for one of the N
%           elliptical contours.
%
%   Written 7/21/2025 by Jess Haley in MATLAB R2024a.
%
%   See also FITGMDIST, EIG.

% Number of contours to draw
numAlpha = length(alpha);

% Chooses points
theta = linspace(0,2*pi,numPoints)';

% Calculates principal directions(PD) and variances (PV)
[PD,PV] = eig(sigma); 

% Construct ellipse fitting 1 standard deviation
ellipseSTD = [cos(theta),sin(theta)]*sqrt(PV)*PD';

% Replicate ellipse for each standard deviation alpha
ellipses = repmat(mu,numPoints,numAlpha) + ...
    repmat(ellipseSTD,1,numAlpha).*repmat(alpha(floor(1:.5:numAlpha+.5)),numPoints,1);

end