function [model, modelMets, solFinalVals]= newCheckMets(solFinalVals,model,biomass_ind)

if exist('biomass_ind', 'var')
    model = removeRxns(model, model.rxns(biomass_ind));
    solFinalVals(biomass_ind)=[];
end
rxnIndsActive = find(abs(solFinalVals)>10^-6);
modelMets.rxnsActive = model.rxns(rxnIndsActive);
modelMets.rxnsInactive = model.rxns(find(abs(solFinalVals)<=10^-6));
sActiveTemp = model.S(:,rxnIndsActive);
metIndsActive = find(logical(sum(abs(sActiveTemp),2)));
metsCarbon = find(~cell2mat(cellfun(@isempty,strfind(model.metFormulas,'C'),'UniformOutput',false)));
modelMets.metIndsActiveCarbon = intersect(metsCarbon,metIndsActive);
modelMets.metsActiveCarbon = model.mets(modelMets.metIndsActiveCarbon);
sActiveCarbonIn = model.S(modelMets.metIndsActiveCarbon,rxnIndsActive);
modelMets.fluxesActive = solFinalVals(rxnIndsActive);
fluxSigns = sign(modelMets.fluxesActive);
modelMets.sActiveDemandCarbDir = bsxfun(@times,sActiveCarbonIn',fluxSigns)';

modelMets.metNames = modelMets.metsActiveCarbon;

% THIS PART OF THE CODE WAS MEANT TO MAKE THE NEW MODEL
% COMPATIBLE WITH OBTAINED RESULTS
% modelMets.metNames = cell(length(modelMets.metsActiveCarbon),1);
% 
% for i = 1:length(modelMets.metNames)
%     curMet = modelMets.metsActiveCarbon{i};
%     curInd = find(strcmp(curMet,model.mets));
%     if ~isempty(curInd)
%         modelMets.metNames{i} = model.mets{curInd};
%     end
% end
