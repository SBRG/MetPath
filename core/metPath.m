function [calcPaths] = metPath(model, modelMets, metsCurCofInorg, cutoffDistance,cutoffFraction)

% it extract the pathways, calculates the weightings, the statistics and
% generates the perturbation scores.
% the output are:
%   resultsTab = a cell array that could be sorted
%   cRes =  vector used to generate resultsTab, useful for other functions

%Defining different S matrices for normal, cofactor, and currency
%metabolites
s = pruneMatrices(model, modelMets, metsCurCofInorg);

%Pathway calculation
calcPaths = calcPathways(model,modelMets, metsCurCofInorg,s,cutoffDistance ,cutoffFraction);

