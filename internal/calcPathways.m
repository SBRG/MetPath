function [calcPaths,srttime,fnstime] = calcPathways(model,modelMets, metsCurCofInorg,s,cutoffDistance ,cutoffFraction)

% For each metabolite extract the reactions involved in the production and 
% its degradation, estimate their weightings and their levels. The levels are
% meant as the distance from the reaction directly involved in the production
% of the metabolite.

% the user can set the distance of the rxns to take in consideration during
% the pathways extraction (cutoffDistance) and a threshold to filter out that
% reactions that slightly partecipates to the pathway to which it belongs

srttime = clock;
srttime = [srttime(4) srttime(5)] ;
display(strcat('started at ',' ',num2str(srttime(1)),':',num2str(srttime(2))))

if ~exist('cutoffDistance', 'var')
    cutoffDistance = 3;
end

% since longer distance will result in smaller
% fraction and further the removal of corresponding pathway

if ~exist('cutoffFraction', 'var')
    cutoffFraction = 0.05;
end

pathwaysProd = zeros(length(modelMets.metIndsActive),length(model.rxns));
pathwaysDeg = zeros(length(modelMets.metIndsActive),length(model.rxns));
probProd = zeros(length(modelMets.metIndsActive),1);
probDeg = zeros(length(modelMets.metIndsActive),1);
% levelsProd = cell(length(modelMets.metIndsActive),1);
% levelsDeg = cell(length(modelMets.metIndsActive),1);

wb = waitbar(0,'init'); 

for i = 1:length(modelMets.metIndsActive)
    
    waitbar(i / length(modelMets.metIndsActive),wb, ['Pathways extracted: ', num2str(i), '/',num2str(length(modelMets.metIndsActive))]) 
    curMet = modelMets.metsActive{i};

    % Reverse direction, production distance

    if ~isempty(find(strcmp(curMet,metsCurCofInorg.currencyPairsComp(:)), 1))
        % It's currency, use the S with inorganic removed
        curDist = distanceDirectional(s.sNoIno, i, cutoffDistance, 0);
    elseif ~isempty(find(strcmp(curMet,metsCurCofInorg.cofactorPairsComp(:)), 1))
        % It's a cofactor, use the S with currency and inorganic mets removed
        curDist = distanceDirectional(s.sNoCurNoIno, i, cutoffDistance, 0);
    elseif ~isempty(find(strcmp(curMet,metsCurCofInorg.inorganicMetsComp(:)), 1))
        % It's a inorganic met, use the full S
        curDist = distanceDirectional(modelMets.sActiveDemandDir, i, cutoffDistance, 0);
    else
        % It's something else, use S with currency, cofactors and  inorganic mets removed
        curDist = distanceDirectional(s.sNoCurNoCofNoIno, i, cutoffDistance, 0);
    end
    curRxnIndsActive = find(curDist>-1);
    if length(curRxnIndsActive)~=length(curDist)
        rxnsToRemove = union(modelMets.rxnsActive(curDist==-1),modelMets.rxnsInactive);
        modelRed = removeRxns(model,rxnsToRemove,false,false);
    else
        rxnsToRemove = modelMets.rxnsInactive;
        modelRed = removeRxns(model,rxnsToRemove,false,false);
    end
    
    
    % Calculate the reduced model imbalances
    curFluxes = modelMets.fluxesActive(curRxnIndsActive);
    curBalances = modelRed.S*curFluxes;
    tol = 10^-8;
    curBalances(abs(curBalances)<tol) = 0;
    curMetIndsUnbalanced = find(curBalances);
    for j = 1:length(curMetIndsUnbalanced)
        curMetDemand = modelRed.mets{curMetIndsUnbalanced(j)};
        modelRed = addReactionNoDup(modelRed,['DM_' curMetDemand],{curMetDemand},1,1,-1000,1000,0,{},{},{},{},false);
        curFluxes = [curFluxes;-curBalances(curMetIndsUnbalanced(j))];

    end
    % Check that the model and flux are now balanced
    max(abs(modelRed.S*curFluxes));

    % Calculate elementary modes
    try
        [P,w] = em_decomp(curFluxes,modelRed);
        % Use weightings to filter out tiny values
        curDemandName = ['DM_' curMet];
        % finding where the P matrix has no demand reaction
        rxnNameInd = strcmp(curDemandName,modelRed.rxns);
        pDemand = find(P(rxnNameInd,:)~=0);
        pWeighted = bsxfun(@times,P(:,pDemand)',w(pDemand)')';
        pSummed = sum(pWeighted,2);
        % Remove components of pSummed that are not at least "cutoffFraction" of the total
        pSummed(abs((pSummed/sum(abs(pSummed))))<cutoffFraction) = 0;
        % Now map this back to the dimensions of model
        for j = 1:length(modelRed.rxns)
            curRxnIndMap = find(strcmp(modelRed.rxns{j},model.rxns));
            if ~isempty(curRxnIndMap)
                pathwaysProd(i,curRxnIndMap) = pSummed(j);
            end
        end
    catch
        probProd(i) = 1;
        for j = 1:length(curRxnIndsActive)
            curRxnIndMap = find(strcmp(modelMets.rxnsActive{curRxnIndsActive(j)},model.rxns));
            if ~isempty(curRxnIndMap)
                pathwaysProd(i,curRxnIndMap) = 1;
            end
        end
    end
    
    
    
%     % adding distance information: distance measurement start from 0
%     
%     ind = match(model.rxns(pathwaysProd(i,:)), modelMets.rxnsActive(curRxnIndsActive));
%     selectedDist = curDist(curDist>-1);
%     selectedDist = selectedDist(ind);
%     
%     if ~isempty(selectedDist)
%         selectedDist = num2str(selectedDist);
%         levelsProd(i) = cell2string(selectedDist);
%     else
%         levelsProd(i) = {''};
%     end
  
    
    
    % Forward direction, degradation distance
    if ~isempty(find(strcmp(curMet,metsCurCofInorg.currencyPairsComp(:)), 1))
        % It's currency, use the full S
        curDist = distanceDirectional(s.sNoIno, i, cutoffDistance, 1);
    elseif ~isempty(find(strcmp(curMet,metsCurCofInorg.cofactorPairsComp(:)), 1))
        % It's a cofactor, use the S with currency and inorganic mets removed
        curDist = distanceDirectional(s.sNoCurNoIno, i, cutoffDistance, 1);
    elseif ~isempty(find(strcmp(curMet,metsCurCofInorg.inorganicMetsComp(:)), 1))
        % It's a inorganic met, use the full S
        curDist = distanceDirectional(modelMets.sActiveDemandDir, i, cutoffDistance, 1);
    else
        % It's something else, use S with currency, cofactors and  inorganic mets removed
        curDist = distanceDirectional(s.sNoCurNoCofNoIno, i, cutoffDistance, 1);
    end

    curRxnIndsActive = find(curDist>-1);
    if length(curRxnIndsActive)~=length(curDist)
        rxnsToRemove = union(modelMets.rxnsActive(curDist==-1),modelMets.rxnsInactive);
        modelRed = removeRxns(model,rxnsToRemove,false,false);
    else
        rxnsToRemove = modelMets.rxnsInactive;
        modelRed = removeRxns(model,rxnsToRemove,false,false);
    end
    
    
    % Calculate the reduced model imbalances
    curFluxes = modelMets.fluxesActive(curRxnIndsActive);
    curBalances = modelRed.S*curFluxes;
    tol = 10^-8;
    curBalances(abs(curBalances)<tol) = 0;
    curMetIndsUnbalanced = find(curBalances);
    for j = 1:length(curMetIndsUnbalanced)
        curMetDemand = modelRed.mets{curMetIndsUnbalanced(j)};
        modelRed = addReactionNoDup(modelRed,['DM_' curMetDemand],{curMetDemand},1,1,-1000,1000,0,{},{},{},{},false);
        curFluxes = [curFluxes;-curBalances(curMetIndsUnbalanced(j))];
    end

    % Calculate elementary modes
    try
        [P,w] = em_decomp(curFluxes,modelRed);
        % Use weightings to filter out tiny values        
        curDemandName = ['DM_' curMet];
        % finding where the P matrix has no demand reaction
        rxnNameInd = strcmp(curDemandName,modelRed.rxns);
        pDemand = find(P(rxnNameInd,:)~=0);
        pWeighted = bsxfun(@times,P(:,pDemand)',w(pDemand)')';
        pSummed = sum(pWeighted,2);
        % Remove components of pSummed that are not at least 5% of the total
        pSummed(abs((pSummed/sum(abs(pSummed))))<cutoffFraction) = 0;
        % Now map this back to the dimensions of model
        for j = 1:length(modelRed.rxns)
            curRxnIndMap = find(strcmp(modelRed.rxns{j},model.rxns));
            if ~isempty(curRxnIndMap)
                pathwaysDeg(i,curRxnIndMap) = pSummed(j);
            end
        end
    catch
        probDeg(i) = 1;
        for j = 1:length(curRxnIndsActive)
            curRxnIndMap = find(strcmp(modelMets.rxnsActive{curRxnIndsActive(j)},model.rxns));
            if ~isempty(curRxnIndMap)
                pathwaysDeg(i,curRxnIndMap) = 1;
            end
        end
    end
    
%     % adding distance information: distance measurement start from 0
% 
%     ind = match(model.rxns(pathwaysDeg(i,:)), modelMets.rxnsActive(curRxnIndsActive));
%     selectedDist = curDist(curDist>-1);
%     selectedDist = selectedDist(ind);
%     
%     if ~isempty(selectedDist)
%         selectedDist = num2str(selectedDist);
%         levelsDeg(i) = cell2string(selectedDist);
%     else
%         levelsDeg(i) = {''};
%     end

   
end

close(wb)
fnstime = clock;
fnstime = [fnstime(4) fnstime(5)] ;
display(strcat('ended at: ',num2str(fnstime(1)),':',num2str(fnstime(2))))

% calcPaths.levelsDeg = levelsDeg;
% calcPaths.levelsProd = levelsProd;
calcPaths.pathwaysDeg = pathwaysDeg;
calcPaths.pathwaysProd = pathwaysProd;

