function [model,modelMets] = defineMets(model, biomass_ind, FBAmode, currencyPairs, cofactorPairs, compartments, allowLoops, inorganicMets)

% this function convert the model in order to be handled by the toolbox,
% define the flux state and generate modelMets struct object necessary for
% the next steps.


% FBAmode permits to set how to solve the LP problem: 1-pFBA, 2-low tolerance,
% 3-minimize fluxes fixing the lb of obj, 4- minimize fluxes

% allowLoops: 1 = allow or 0 = not allow

% inorganicMets: define the list of inorganic metabolites to take in
% consideration for pathways extraction. If not present the default list
% will be used

model = convertModel(model);

try
    model = calculateFluxState(model,FBAmode,allowLoops);
catch
    model = calculateFluxState(model);
end

try
    [model,modelMets] = checkMets(model,biomass_ind, inorganicMets, currencyPairs, cofactorPairs, compartments);
catch
    [model,modelMets] = checkMets(model, biomass_ind, 0, currencyPairs, cofactorPairs, compartments);
end