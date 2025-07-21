function [estimatedPosterior] = estimatePosterior(model,features,X)
%[estimatedPosterior] = ESTIMATEPOSTERIOR(model,features,X)
%
%   ESTIMATEPOSTERIOR estimates the probability of observations of a 
%   feature X belonging to the first cluster in a discriminant analysis
%   model by marginalizing the conditional probabilities over the remaining
%   features Y and Z and numerically integrating using an adaptive 
%   quadrature method over the product of the QDA estimated conditional 
%   probability distribution and the kernel density estimate of the joint 
%   probability distribution.
%
%   INPUTS:
%       - model [obj]: 2-cluster discriminant analysis classifier object;
%           created with LABELENCOUNTERS
%       - features [double]: matrix containing all observations of three 
%           data variables (i.e. [X, Y, Z])
%       - X [double]: array containing observations with only the first
%           variable
%
%   OUTPUTS:
%       - estimatedPosterior [double]: posterior probability of
%           observations of X belonging to cluster 1
%
%   Written 3/26/2024 by Jess Haley in MATLAB R2024a.
%
%   See also LABELENCOUNTERS, INTEGRAL, MVKSDENSITY, PREDICT.

% Check that features: rows=observations, cols=variables
if ndims(features) == 2 & diff(size(features)) > 0
    features = features';
end

% P(cluster1|X) = P(cluster1,X) / P(X)
estimatedPosterior = integrateProb_dY(model,features,X)./kdeN(features(:,1),X);

    % P(cluster1,X) = integral P(cluster1|X,Y) * P(X,Y) dY
    function [pC1_and_X] = integrateProb_dY(model,features,X)
        pC1_and_X = integral(@(Y) integrateProb_dZ(model,features,X,Y),...
            -Inf,Inf,'ArrayValue',true,'AbsTol',1e-4,'RelTol',1e-4);
    end

    % P(cluster1,X,Y) = integral P(cluster1|X,Y,Z) * P(X,Y,Z) dZ
    function [pC1_and_XY] = integrateProb_dZ(model,features,X,Y)
        pC1_and_XY = integral(@(Z) clusterProb(model,X,Y,Z).*observationProb(features,X,Y,Z),...
            -Inf,Inf,'ArrayValue',true,'AbsTol',1e-4,'RelTol',1e-4);
    end

    % P(cluster1|X,Y,Z)
    function [pC1_given_XYZ] = clusterProb(model,X,Y,Z)
        [~,posteriorProb] = predict(model,[X,Y,Z]);
        pC1_given_XYZ = posteriorProb(:,1);
    end

    % P(X,Y,Z)
    function [pXYZ] = observationProb(features,X,Y,Z)
        pXYZ = kdeN(features,[X,Y,Z]); % P(X,Y,Z)
    end

    % N-Dimensional Kernel Density Estimate
    function [pdf] = kdeN(features,points)
        d = size(features,2); % # dimensions
        n = size(features,1); % # observations
        bw = std(features,[],1,'omitnan').*(4/(n*(d+2))).^(1/(d+4)); % Silverman's rule
        pdf = mvksdensity(features,points,'Bandwidth',bw);
    end 

end