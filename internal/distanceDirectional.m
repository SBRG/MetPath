function distVec = distanceDirectional(sDir, metInd, depth,  forwardFlag)

%This function starts at one particular metabolite, and given a direction
%stoichiometric matrix, being one that has been sign corrected to be
%irreversible with respect to a particular flux distribution, finds the
%distance forward and reverse between a hypothetical demand reaction on the
%metabolite and other reactions in the network, to a pre-asssigned depth.

%Need to delete the reactions that degrade the metabolite for the
%production distance, and delete the reactions that produce the metabolite
%for the degradation distance, to prevent cycling

rxnIndsCons = find(sDir(metInd,:)<0);
rxnIndsProd = find(sDir(metInd,:)>0);

sPosDeg = sDir;
sPosDeg(sPosDeg<0) = 0;
sDirDeg = sDir;
sDirDeg(:,rxnIndsProd) = 0;

sNegProd = sDir;
sNegProd(sNegProd>0) = 0;
sDirProd = sDir;
sDirProd(:,rxnIndsCons) = 0;

distVec = -1*ones(size(sDir,2),1);

if forwardFlag == 1
    %Start with the reactions that consume the metabolite of interest
    distVec(rxnIndsCons) = 0;
    curRxnInds = rxnIndsCons;
    for i = 1:depth
        %Get metabolites that are produced by reactions at the current
        %level
        curMetInds = find(logical(sum(abs(sPosDeg(:,curRxnInds)),2)));
        %Get reactions that consume these metabolites
        sNew = sDirDeg(curMetInds,:);
        sNew(sNew>0) = 0;
        curRxnInds = find(logical(sum(abs(sNew),1)));
        for j = 1:length(curRxnInds)
            if distVec(curRxnInds(j))==-1
                distVec(curRxnInds(j)) = i;
            end
        end
    end
else
    %Start with the reactions that produce the metabolite of interest
    distVec(rxnIndsProd) = 0;
    curRxnInds = rxnIndsProd;
    for i = 1:depth
        %Get metabolites that are consumed by reactions at the current
        %level
        curMetInds = find(logical(sum(abs(sNegProd(:,curRxnInds)),2)));
        %Get reactions that produce these metabolites
        sNew = sDirProd(curMetInds,:);
        sNew(sNew<0) = 0;
        curRxnInds = find(logical(sum(abs(sNew),1)));
        for j = 1:length(curRxnInds)
            if distVec(curRxnInds(j))==-1
                distVec(curRxnInds(j)) = i;
            end
        end
    end
end

