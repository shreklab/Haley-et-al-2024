function [significance,pAdjusted] = benjaminiHochberg(pValues,Q)
% [SIGNIFICANCE,PADJUSTED] = BENJAMINIHOCHBERG(PVALUES,Q)
%
%   BENJAMINIHOCHBERG corrects for multiple comparisons by controlling the
%   false discovery rate (FDR) using the Benjamini-Hochberg procedure. The
%   function identifies which p-values in a set are significant given an
%   FDR, Q, and calculates their corresponding adjusted p-values.
%
%   INPUTS:
%       - pValues [double]: An array of p-values from multiple hypothesis
%           tests.
%       - Q [double]: The acceptable false discovery rate (e.g., 0.05).
%
%   OUTPUTS:
%       - significance [logical]: A logical array of the same size as
%           pValues where 'true' indicates a significant result after
%           FDR correction.
%       - pAdjusted [double]: An array of the adjusted p-values, which can
%           be directly compared against Q to determine significance.
%
%   Written 7/11/2025 by Jess Haley in MATLAB R2024a.
%
%   See also SORT, RESHAPE, FILLMISSING.

pValues = reshape(pValues,1,[]);

% Sort p-values (ascending)
[pValuesSorted,indSort] = sort(pValues,'ascend');

% Compute Benjamini-Hochbery critical value (i/m)Q
m = length(pValues);
i = 1:m;
imQ = Q*i./m;

% Get adjusted p-values
pAdjustedSorted = m*pValuesSorted./i;

% Get index of large p-value smaller than it's BH equivalent (both
% statements are equivalent)
maxBH = find(pValuesSorted(:) < imQ(:),1,'last');
% maxBH =  find(pAdjustedSorted(:) < Q,1,'last');

% Get signficance w/ BH-correction
significance = false(size(pValues));
significance(indSort(1:maxBH)) = true;

% Ensure adjusted p-values are all ascending (replace with next lowest
% adjusted p-value if descending)
pAdjustedSorted(diff([0,pAdjustedSorted]) < 0) = NaN;
pAdjustedSorted = fillmissing(pAdjustedSorted,'previous');
[~,~,indUnsort] = intersect(1:m,indSort);
pAdjusted = pAdjustedSorted(indUnsort);

end