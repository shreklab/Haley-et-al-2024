function [imageFit,parameters] = getFluorescenceBackground(image,method,exposureTime,startPoint)
% [imageFit,parameters] = GETFLUORESCENCEBACKGROUND(image,method,exposureTime,startPoint)
%
%   GETFLUORESCENCEBACKGROUND takes in an image of an empty agar plate and
%   finds a filtered/fitted landscape which highlights the diffusion of the
%   light source across the field of views. The landscape can be fit by
%   assuming the light diffuses as an elliptic paraboloid OR by smoothing
%   the image via an average filtering operation.
%
%   INPUTS:
%       - image [MxN num]: image of an empty agar plate
%       - method [str]: fitting method to use (i.e. 'ellipse' OR 'smooth')
%       - exposureTime [double]: exposure time of the image in ms
%       - startPoint [1x5 double]: start point of each of the coefficients
%           [a,b,x0,y0,z0] in the elliptic paraboloid model
%
%   OUTPUTS:
%       - imageFit [MxN double]: image of the filtered/fit light source
%       - parameters [table]: table of parameters
%           - rsquare [double]: r-squared value
%           - numOutliers [double]: number of outliers replaced ('smooth'
%               only)
%           - a [double]: paraboloid coefficient for x ('ellipse' only)
%           - b [double]: paraboloid coefficient for y ('ellipse' only)
%           - x0 [double]: paraboloid origin for x ('ellipse' only)
%           - y0 [double]: paraboloid origin for y ('ellipse' only)
%           - z0 [double]: paraboloid origin for z ('ellipse' only)
%    
%   Written 1/29/2024 by Jess Haley in MATLAB R2023b.
%
%   See also ANALYZEGFP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Fit elliptic paraboloid to background image

if strcmp(method,'ellipse')

    % Define elliptic paraboloid model (z = a*(x-x0)^2 + b*(y-y0)^2 + z0)
    ellipticParaboloid = fittype('a*(x-x0)^2 + b*(y-y0)^2 + z0',...
        dependent='z',independent={'x','y'},...
        coefficients={'a','b','x0','y0','z0'});

    % Get x, y, and z for model
    [x,y] = find(ones(size(image)));
    z = double(image);

    % Set exposure time (if undefined)
    if isempty(exposureTime)
        exposureTime = 1000;
    end

    % Set model start point (if undefined)
    if isempty(startPoint)
        startPoint = [-3e-7*exposureTime,-3e-7*exposureTime,...
            max(x)/2,max(y)/2,prctile(z,99,'all')];
    end

    % Fit elliptic paraboloid to image
    [fitObject,gof,~] = fit([x,y],z(:),ellipticParaboloid,...
        'StartPoint',startPoint);

    % Get image from fit object
    imageFit = reshape(fitObject(x,y),size(z));

    % Save parameters
    parameters.a = fitObject.a;
    parameters.b = fitObject.b;
    parameters.x0 = fitObject.x0;
    parameters.y0 = fitObject.y0;
    parameters.z0 = fitObject.z0;
    parameters.rsquare = gof.rsquare;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Moving average of background image

elseif strcmp(method,'smooth')

    % Set image data to double
    z = double(image);

    % Find and replace outliers
    outliers = isoutlier(z,'movmedian',400,1) & ...
        isoutlier(z,'movmedian',400,2);
    zClean = z;
    zClean(outliers) = NaN;
    zClean = fillmissing(zClean,'linear');

    % Smooth image
    imageFit = imfilter(zClean,fspecial('average',20),'symmetric');

    % Save parameters
    parameters.numOutliers = sum(outliers,'all');
    parameters.rsquare = 1 - sum((z - imageFit).^2)/sum((z - mean(z,'all')).^2);

end

parameters = struct2table(parameters);

end