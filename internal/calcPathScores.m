function scores = calcPathScores(model, modelMets, fMap, calcPaths, numPerms)

% score the production scores, degradation scores and the aggregate 
% perturbation score using the extracted pathways

%Pathways for production/degradation
pProdString = cell(length(modelMets.metIndsActive),1);
pDegString = cell(length(modelMets.metIndsActive),1);

%Weightings for production/degradation
wProdString = cell(length(modelMets.metIndsActive),1);
wDegString = cell(length(modelMets.metIndsActive),1);

%Differential gene expression score for the pathway
scoresProd = NaN*ones(length(modelMets.metIndsActive),1);
scoresDeg = NaN*ones(length(modelMets.metIndsActive),1);

%Permutation p-value for the pathway being down-regulated
pValLowProd = NaN*ones(length(modelMets.metIndsActive),1);
pValLowDeg = NaN*ones(length(modelMets.metIndsActive),1);

%Permutation p-value for the pathway being up-regulated
pValHighProd = NaN*ones(length(modelMets.metIndsActive),1);
pValHighDeg = NaN*ones(length(modelMets.metIndsActive),1);

for i = 1:length(modelMets.metIndsActive)
    %% prod
    curPathProd = calcPaths.pathwaysProd(i,:);
    curRxnsActive = model.rxns(find(curPathProd));
    curRxnsMap = intersect(curRxnsActive, fMap.corrRxnUniqueData);
    curRxnIndsMap = cell2mat(cellfun(@(x) find(strcmp(x,model.rxns)),curRxnsMap,'UniformOutput',false));
    curWeights = abs(curPathProd(curRxnIndsMap)/sum(abs(curPathProd(curRxnIndsMap))));
    
    curCell = curRxnsMap;
    curString = '';
    if ~isempty(curCell)
        for j = 1:(length(curCell)-1)
            curString = [curString num2str(curCell{j}) ';'];
        end
        curString = [curString num2str(curCell{end})];
    end
    pProdString{i} = curString;
    
    curMat = curWeights;
    curString = '';
    if ~isempty(curCell)
        for j = 1:(length(curCell)-1)
            curString = [curString num2str(curMat(j)) ';'];
        end
        curString = [curString num2str(curMat(end))];
    end
    wProdString{i} = curString;
    
    curFoldInds = cell2mat(cellfun(@(x) find(strcmp(x,fMap.corrRxnUniqueData)),curRxnsMap,'UniformOutput',false));
    if ~isempty(curFoldInds)
        curFolds = fMap.foldChangeRxnData(curFoldInds);
        curScoreProd = nansum(curWeights.*curFolds'); % replaced sum with nansum
        scoresProd(i) = curScoreProd;
        [pValLow, pValHigh] = pathStats(curScoreProd, curWeights, fMap.foldChangeRxnData, numPerms);
        pValLowProd(i) = pValLow;
        pValHighProd(i) = pValHigh;
    end
    
    % add subsystem info
    if ~isempty(curCell)
        tmp_curSS = model.subSystems(match(curCell, model.rxns));
        if ~isempty(tmp_curSS)
            pSS(i,1) = cell2string(tmp_curSS);
        end
    end
    

    %% deg
    curPathDeg = calcPaths.pathwaysDeg(i,:);
    curRxnsActive = model.rxns(find(curPathDeg));
    curRxnsMap = intersect(curRxnsActive, fMap.corrRxnUniqueData);
    curRxnIndsMap = cell2mat(cellfun(@(x) find(strcmp(x,model.rxns)),curRxnsMap,'UniformOutput',false));
    curWeights = abs(curPathDeg(curRxnIndsMap)/sum(abs(curPathDeg(curRxnIndsMap))));
    curFoldInds = cell2mat(cellfun(@(x) find(strcmp(x,fMap.corrRxnUniqueData)),curRxnsMap,'UniformOutput',false));
    if ~isempty(curFoldInds)
        curFolds = fMap.foldChangeRxnData(curFoldInds);
        curScoreDeg = nansum(curWeights.*curFolds'); % replaced sum with nansum,
        scoresDeg(i) = curScoreDeg;
        [pValLow, pValHigh] = pathStats(curScoreDeg, curWeights, fMap.foldChangeRxnData, numPerms);
        pValLowDeg(i) = pValLow;
        pValHighDeg(i) = pValHigh;
    end
    
    curCell = curRxnsMap;
    curString = '';
    if ~isempty(curCell)
        for j = 1:(length(curCell)-1)
            curString = [curString num2str(curCell{j}) ';'];
        end
        curString = [curString num2str(curCell{end})];
    end
    pDegString{i} = curString;
    
    curMat = curWeights;
    curString = '';
    if ~isempty(curCell)
        for j = 1:(length(curCell)-1)
            curString = [curString num2str(curMat(j)) ';'];
        end
        curString = [curString num2str(curMat(end))];
    end
    wDegString{i} = curString;
    
    % subsystem info
    if ~isempty(curCell)
        tmp_curSS = model.subSystems(match(curCell, model.rxns));
        if ~isempty(tmp_curSS)
            dSS(i,1) = cell2string(tmp_curSS);  
        end
    end
    
end


scores.pProdString = pProdString;
scores.wProdString = wProdString;
scores.scoresProd = scoresProd;
scores.pDegString = pDegString;
scores.wDegString = wDegString;
scores.scoresDeg = scoresDeg;
scores.dSubSyst=dSS;
scores.pSubSyst=pSS;
scores.pValLowProd = pValLowProd;
scores.pValHighProd = pValHighProd;
scores.pValLowDeg = pValLowDeg;
scores.pValHighDeg = pValHighDeg;
