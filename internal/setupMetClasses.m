function [metsCurCofInorg] = setupMetClasses(currencyPairs, cofactorPairs, compartments, inorganicMets)

% This function
% Inputs
%   currencyPairs: Pairs of currency metabolites to ignore in pathway calc
%   cofactorPairs: Pairs of cofactor metabolits to ignore in pathway calc
%   compartments: Compartments to search for currency/cofactor mets
% Outputs

%% Set up list of currency and cofactor metabolites

%THIS FUNCTION COULD PROBABLY BE SEPARATELY DONE AT THE START OF THE
%WORKFLOW
%Assembled list of cofactors
metsCurCofInorg.allPairs = [currencyPairs;cofactorPairs];

metsCurCofInorg.currencyPairsComp = {};
metsCurCofInorg.currencyPairsMap = [];
for i = 1:length(currencyPairs)
    for j = 1:length(compartments)
        metsCurCofInorg.currencyPairsComp = [metsCurCofInorg.currencyPairsComp;{[currencyPairs{i,1} compartments{j}],[currencyPairs{i,2} compartments{j}]}];
        metsCurCofInorg.currencyPairsMap = [metsCurCofInorg.currencyPairsMap;i];
    end
end

metsCurCofInorg.cofactorPairsComp = {};
metsCurCofInorg.cofactorPairsMap = [];
for i = 1:length(cofactorPairs)
    for j = 1:length(compartments)
        metsCurCofInorg.cofactorPairsComp = [metsCurCofInorg.cofactorPairsComp;{[cofactorPairs{i,1} compartments{j}],[cofactorPairs{i,2} compartments{j}]}];
        metsCurCofInorg.cofactorPairsMap = [metsCurCofInorg.cofactorPairsMap;i];
    end
end


% % adding the currency metabolites
% if ~isempty(currencyPairs)
%     metsCurCofInorg.currency=[currencyPairs(:,1);currencyPairs(:,2)];
% else 
%     metsCurCofInorg.currency = {};
% end
% 
% %SHOULDN'T HAVE TO HARD CODE CO2 IN HERE
% % metsCurCofInorg.currency=[metsCurCofInorg.currency;{'co2'}];
% metsCurCofInorg.currencyComp = {};
% metsCurCofInorg.currencyMap = [];
% for i = 1:length(metsCurCofInorg.currency)
%     for j = 1:length(compartments)
%         metsCurCofInorg.currencyComp = [metsCurCofInorg.currencyComp;{[metsCurCofInorg.currency{i,1} compartments{j}]}];
%         metsCurCofInorg.currencyMap = [metsCurCofInorg.currencyMap;i];
%     end
% end

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

metsCurCofInorg.inorganicMetsComp = inorganicMetsComp;
metsCurCofInorg.inorganicMetsMap = inorganicMetsMap;