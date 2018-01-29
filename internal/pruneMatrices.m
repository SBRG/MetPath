function s = pruneMatrices(model, modelMets)

% generate different S matrix for each type of metabolite

s.sNoCurNoIno = modelMets.sActiveDemandDir;
s.sNoCurNoCofNoIno = modelMets.sActiveDemandDir;
s.sNoIno = modelMets.sActiveDemandDir;


for i = 1:length(modelMets.rxnsActive)
    
    curRxn = modelMets.rxnsActive{i};
    curRxnInd = find(strcmp(curRxn,model.rxns));
    curRxnActiveInd = find(strcmp(curRxn,modelMets.rxnsActive));
    curMets = model.mets(find(model.S(:,curRxnInd)));
    curMetFormulas = model.metFormulas(find(model.S(:,curRxnInd)));
    numCarb = length(find(findregexp(curMetFormulas, 'C\w[^a-z]',1)));
    
    numCur = 0;
    curCurMets = {};
    for j = 1:length(model.curCof.currencyPairsMap)
        if ~isempty(find(strcmp(model.curCof.currencyPairsComp{j,1},curMets))) && ~isempty(find(strcmp(model.curCof.currencyPairsComp{j,2},curMets)))
            numCur = numCur + 2;
            curCurMets = [curCurMets;model.curCof.currencyPairsComp{j,1};model.curCof.currencyPairsComp{j,2}];
        end
    end


    numCof = 0;
    model.curCofMets = {};
    for j = 1:length(model.curCof.cofactorPairsMap)
        if ~isempty(find(strcmp(model.curCof.cofactorPairsComp{j,1},curMets))) && ~isempty(find(strcmp(model.curCof.cofactorPairsComp{j,2},curMets)))
            numCof = numCof + 2;
            model.curCofMets = [model.curCofMets;model.curCof.cofactorPairsComp{j,1};model.curCof.cofactorPairsComp{j,2}];
        end
    end

    
    % different metabolites type use different matrices:
    % currencies -> no ino
    % inorganicMets -> full S
    % cofactors -> no cur no ino
    % other mets -> no cur no cof no ino
    
    % Generate another list with special Mets to prune noCur->noCurNoIno and
    % noCurNoCof-> noCurNoCofNoIno

    % list of inorganic mets
    curInoMets = {};
    for j = 1:length(model.curCof.inorganicMetsComp)
        if ~isempty(find(strcmp(model.curCof.inorganicMetsComp{j,1},curMets)))
            curInoMets = [curInoMets;model.curCof.inorganicMetsComp{j,1}];
        end
    end    

    
    % Pruning Currencies
    % If there's currency and more carbon than currency, remove currency
    % from the reaction in the noCurNoIno and noCurNoCofNoIno
    % matrices
    remanining_c1 = numCarb;
    if ~isempty(curCurMets) && (numCarb-numCur) > 0
        remanining_c1 = remanining_c1-numCur;
        for j = 1:length(curCurMets)
            curInd = find(strcmp(curCurMets{j},modelMets.metsActive));
            s.sNoCurNoIno(curInd,curRxnActiveInd) = 0; 
            s.sNoCurNoCofNoIno(curInd,curRxnActiveInd) = 0;
            
        end
    end
    
    % Pruning inorganic Mets
	% If there's a inorganic mets, remove them from the rxns in the 
    % noCurNoIno, noIno and noCurNoCofNoIno
    if ~isempty(curInoMets)
        for j = 1:length(curInoMets)
            curInd = find(strcmp(curInoMets{j},modelMets.metsActive));
            s.sNoCurNoIno(curInd,curRxnActiveInd) = 0; 
            s.sNoCurNoCofNoIno(curInd,curRxnActiveInd) = 0;
            s.sNoIno(curInd,curRxnActiveInd) = 0;
        end
    end

        
    % Pruning Cofactors
    % If there's currency or cofactor and remaining carbon, remove
    % currency/cofactors from the reaction in the NoCurNoCof S matrix
    nummodel.curCof = numCur+numCof;
    if ~isempty(model.curCofMets) && (numCarb-nummodel.curCof) > 0
        for j = 1:length(model.curCofMets)
                curInd = find(strcmp(model.curCofMets{j},modelMets.metsActive));
                s.sNoCurNoCofNoIno(curInd,curRxnActiveInd) = 0;
        end
    end
    
    
    % removing single currency
    % finding currencies:
    numCur_single=0;
    curCurMets_single={};
    for j = 1:length(model.curCof.currencyMap)
        if ~isempty(find(strcmp(model.curCof.currencyComp{j,1},curMets)))
            numCur_single = numCur_single + 1;
            curCurMets_single = [curCurMets_single;model.curCof.currencyComp{j,1}];
        end
    end
    
    % If there's currency and still remaining carbon remove it from the reaction
    % in both S matrices
    if ~isempty(curCurMets_single) && (remanining_c1-numCur_single) > 0
        for j = 1:length(curCurMets_single)
            curInd = find(strcmp(curCurMets_single{j},modelMets.metsActive));
            s.sNoCurNoIno(curInd,curRxnActiveInd) = 0; 
            s.sNoCurNoCofNoIno(curInd,curRxnActiveInd) = 0;
        end
    end
    
    
end
