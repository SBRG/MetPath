function [modelNoBM, structActive, nonCarbonMets, fluxesRed]= getActiveNetwork(model,biomassInd, fluxes, inorganicMets, compartments,tolFlux)

% This function takes a flux state and returns the active sets of
% metabolites and reactions
% Inputs
%   model: Constraint-based model to be used
%   biomassInd: Index of the biomass objective reaction
%   fluxes: Representative flux solution
%   inorganicMets: 
%      1 for a default list:
%       {'o2', 'so2','so3','so4','nh4','no2','no3','fe2','fe3', 'h2o',...
%        'co2','co','h2o2','o2s','h2s','etha', 'no','fe3hox','3fe4s',...
%        '4fe4s','2fe2s', 'etoh','mobd','cu','cu2'};
%      0 for an empty list
%      Or a user-provided list e.g. {'no','fe2'}
% Outputs
%   modelNoBM: Model with biomass reactions removed to prevent non-pathway
%   connections between metabolites
%   structActive: Structure of active and inactive metabolites and reactions
%   nonCarbonMets: Metabolites not containing any carbon
%   fluxesRed: Corresponding flux solution with biomass removed

%% Setup undesirable fluxes and reactions

% inorganic S
%Allow two default options
if ~iscell(inorganicMets)
    if inorganicMets == 1
        inorganicMets = {'o2', 'so2','so3','so4','nh4','no2','no3','fe2','fe3', 'h2o','co2','co','h2o2','o2s','h2s',...
            'etha', 'no','fe3hox','3fe4s','4fe4s','2fe2s', 'etoh','mobd','cu','cu2'};
    elseif inorganicMets == 0
        inorganicMets = {};
    end
end

inorganicMetsMap = [];
inorganicMetsComp= {};
for i = 1:length(inorganicMets)
    for j = 1:length(compartments)
        inorganicMetsComp = [inorganicMetsComp;{[inorganicMets{i} compartments{j}]}];
        inorganicMetsMap = [inorganicMetsMap;i];
    end
end

%All non-carbon metabolites
nonCarbonMets = model.mets(find(~findregexp(model.metFormulas, 'C\w[^a-z]',1)));

% Remove the biomass reaction for pathway calculation
fluxesRed = fluxes;
modelNoBM = model;
if exist('biomassInd', 'var')
    modelNoBM = removeRxns(model, model.rxns(biomassInd));
    fluxesRed(biomassInd)=[];
end

%% Check which metabolites are part of flux-carrying reactions

%Flux-carrying reactions and S matrix
rxnIndsActive = find(abs(fluxesRed)>tolFlux);
structActive.rxnsActive = modelNoBM.rxns(rxnIndsActive);
structActive.rxnsInactive = modelNoBM.rxns(find(abs(fluxesRed)<=tolFlux));
sActiveTemp = modelNoBM.S(:,rxnIndsActive);

%Metabolites in flux-carrying reactions
metIndsActive = find(logical(sum(abs(sActiveTemp),2)));

%Restrict to carbon-carrying metabolites
metsCarbon = find(findregexp(modelNoBM.metFormulas, 'C\w[^a-z]',1));

%Add in any additional inorganic mets to be used
usedMets = [metsCarbon; find(ismember(modelNoBM.mets, inorganicMetsComp))];

%Find subset of active metabolites that are also carbon containing or on
%the intended inorganic list
structActive.metIndsActive = intersect(usedMets,metIndsActive);
structActive.metsActive= modelNoBM.mets(structActive.metIndsActive);
structActive.metNames = structActive.metsActive;

%Resulting S matrix and fluxes from reaction and metabolite restrictions
sActiveIn = modelNoBM.S(structActive.metIndsActive,rxnIndsActive);
structActive.fluxesActive = fluxesRed(rxnIndsActive);
fluxSigns = sign(structActive.fluxesActive);

%Directional S matrix such that all reactions are proceeding in the forward
%direction
structActive.sActiveDemandDir = bsxfun(@times,sActiveIn',fluxSigns)';

