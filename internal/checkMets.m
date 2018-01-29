function [model, modelMets, others]= checkMets(model,biomass_ind, inorganicMets, currencyPairs, cofactorPairs, compartments)


curCof.allPairs = [currencyPairs;cofactorPairs];

curCof.currencyPairsComp = {};
curCof.currencyPairsMap = [];
for i = 1:length(currencyPairs)
    for j = 1:length(compartments)
        curCof.currencyPairsComp = [curCof.currencyPairsComp;{[currencyPairs{i,1} compartments{j}],[currencyPairs{i,2} compartments{j}]}];
        curCof.currencyPairsMap = [curCof.currencyPairsMap;i];
    end
end

curCof.cofactorPairsComp = {};
curCof.cofactorPairsMap = [];
for i = 1:length(cofactorPairs)
    for j = 1:length(compartments)
        curCof.cofactorPairsComp = [curCof.cofactorPairsComp;{[cofactorPairs{i,1} compartments{j}],[cofactorPairs{i,2} compartments{j}]}];
        curCof.cofactorPairsMap = [curCof.cofactorPairsMap;i];
    end
end


% adding the currency metabolites
curCof.currency=[currencyPairs(:,1);currencyPairs(:,2)];

curCof.currency=[curCof.currency;{'co2'}];
curCof.currencyComp = {};
curCof.currencyMap = [];
for i = 1:length(curCof.currency)
    for j = 1:length(compartments)
        curCof.currencyComp = [curCof.currencyComp;{[curCof.currency{i,1} compartments{j}]}];
        
        curCof.currencyMap = [curCof.currencyMap;i];
    end
end




%% inorganic S

if ~exist('inorganicMets', 'var')
    inorganicMets = {'o2', 'so2','so3','so4','nh4','no2','no3','fe2','fe3', 'h2o','co2','co','h2o2','o2s','h2s',...
        'etha', 'no','fe3hox','3fe4s','4fe4s','2fe2s', 'etoh','mobd','cu','cu2'};
elseif inorganicMets == 0
    inorganicMets = {};
end

others = model.mets(find(~findregexp(model.metFormulas, 'C\w[^a-z]',1)));

inorganicMetsMap = [];
inorganicMetsComp= {};
for i = 1:length(inorganicMets)
    for j = 1:length(compartments)
        inorganicMetsComp = [inorganicMetsComp;{[inorganicMets{i} compartments{j}]}];
        inorganicMetsMap = [inorganicMetsMap;i];
    end
end

curCof.inorganicMetsComp = inorganicMetsComp;
curCof.inorganicMetsMap = inorganicMetsMap;
%%

if exist('biomass_ind', 'var')
    model = removeRxns(model, model.rxns(biomass_ind));
    model.solFinalVals(biomass_ind)=[];
end

rxnIndsActive = find(abs(model.solFinalVals)>10^-6);
modelMets.rxnsActive = model.rxns(rxnIndsActive);
modelMets.rxnsInactive = model.rxns(find(abs(model.solFinalVals)<=10^-6));
sActiveTemp = model.S(:,rxnIndsActive);
metIndsActive = find(logical(sum(abs(sActiveTemp),2)));
metsCarbon = find(findregexp(model.metFormulas, 'C\w[^a-z]',1));
usedMets = [metsCarbon; find(ismember(model.mets, inorganicMetsComp))];
modelMets.metIndsActive = intersect(usedMets,metIndsActive);
modelMets.metsActive= model.mets(modelMets.metIndsActive);
sActiveIn = model.S(modelMets.metIndsActive,rxnIndsActive);
modelMets.fluxesActive = model.solFinalVals(rxnIndsActive);
fluxSigns = sign(modelMets.fluxesActive);
modelMets.sActiveDemandDir = bsxfun(@times,sActiveIn',fluxSigns)';

modelMets.metNames = modelMets.metsActive

model.curCof = curCof;
