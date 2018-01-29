function [resultsTab, cRes] = metPath(model, modelMets, fMap,cutoffDistance,cutoffFraction)

% it extract the pathways, calculates the weightings, the statistics and
% generates the perturbation scores.
% the output are:
%   resultsTab = a cell array that could be sorted
%   cRes =  vector used to generate resultsTab, useful for other functions

s = pruneMatrices(model, modelMets);


try
    [calcPaths] = calcPathways(model,modelMets,s,cutoffDistance ,cutoffFraction);
catch
    [calcPaths] = calcPathways(model,modelMets,s);
end

cRes = calcRes(model, modelMets, fMap, calcPaths);
cRes.pathwaysDeg = calcPaths.pathwaysDeg;
cRes.pathwaysProd = calcPaths.pathwaysProd;
resultsTab = createResultsTab(modelMets, cRes);