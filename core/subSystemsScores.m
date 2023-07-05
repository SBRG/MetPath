function subSystemsPerturbation = subSystemsScores(model, cRes, modelMets,paths1, model2, cRes2, modelMets2,paths2)


% this function permit by means of perturbation score which are the most
% perturbed subSystems in the model. 
    
utilizationTab = calcAggregateScores(modelMets, cRes);

cnt = {};
strCnt = {};
strScore = {};
paths = utilizationTab(:,1);
score = utilizationTab(:,4);
totCnt = {};
totScore = {};

ssTotCnt = count({model.subSystems});
ssTotCnt = ssTotCnt{:};

for i=1:numel(paths)
 
% counting the subSyst.

    curInd = ismember(modelMets.metNames, paths(i));
    curProd = model.rxns(find(paths1.pathwaysProd(curInd,:)));
    curDeg = model.rxns(find(paths1.pathwaysDeg(curInd,:)));
    curRxnsAct = unique([curProd; curDeg]); 
    curRxnsAct = intersect(curRxnsAct, modelMets.rxnsActive);
    curSubSystems = model.subSystems(match(curRxnsAct, model.rxns)); 
    curCnt = count({curSubSystems});
    
%% standardization

    tmp = curCnt{1};
    if ~isempty(tmp)
        check = 1;
        for n = 1:numel(tmp(:,1))
            ind = find(strcmp(tmp(n,1),ssTotCnt(:,1)));
            tmp(n,2)  = num2cell(cell2mat(tmp(n,2))/cell2mat(ssTotCnt(ind,2)));
        end
        curCnt = {tmp};
    else
        check = 0;
        tmp(1,1) = {''};
        tmp(1,2)  = {0};
        curCnt = {tmp};
    end
    
%% scoring

    cnt = [cnt;curCnt];
    curCnt = curCnt{:};
    curCntScore = curCnt;
    curCntScore(:,2) = num2cell(cell2mat(curCntScore(:,2)) * cell2mat(score(i)));
    curCnt(:,2) = num2cell(cell2mat(curCnt(:,2)) * sign(cell2mat(score(i))));

    
    if check == 1
        curCnt = flipud(sortrows(curCnt,2));
        curCntScore = flipud(sortrows(curCntScore,2));
    end
    tmp = {};
    tmpS ={};    
    
    for k = 1:size(curCnt,1)
        tmp(k) = strcat(curCnt(k,1),' -> ' ,num2str(curCnt{k,2}));
        tmpS(k) = strcat(curCntScore(k,1),' -> ' ,num2str(curCntScore{k,2}));
    end
    
    curStr = append(tmp);
    strCnt = [strCnt;curStr];
    
    curStrScore = append(tmpS);
    strScore = [strScore;curStrScore];
    
    totCnt = [totCnt;curCnt];
    totScore = [totScore;curCntScore];
   
end

topCount = flipud(sortrows(collapseSum(totCnt),2));
topScore = flipud(sortrows(collapseSum(totScore),2));


results1 = topScore;


%% comparing agains second growing condition

if exist('model2', 'var')

    model = model2;
    utilizationTab = calcAggregateScores(modelMets, cRes2);
    modelMets = modelMets2;
    cnt = {};
    strCnt = {};
    strScore = {};
    paths = utilizationTab(:,1);
    score = utilizationTab(:,4);
    totCnt = {};
    totScore = {};

    ssTotCnt = count({model.subSystems});
    ssTotCnt = ssTotCnt{:};

    for i=1:numel(paths)
 
    %   counting the subSyst.

        curInd = ismember(modelMets.metNames, paths(i));
        curProd = model.rxns(find(paths2.pathwaysProd(curInd,:)));
        curDeg = model.rxns(find(paths2.pathwaysDeg(curInd,:)));
        curRxnsAct = unique([curProd; curDeg]); 
        curRxnsAct = intersect(curRxnsAct, modelMets.rxnsActive);
        curSubSystems = model.subSystems(match(curRxnsAct, model.rxns)); 
        curCnt = count({curSubSystems});
    
    %% standardization

    tmp = curCnt{1};
        if ~isempty(tmp)
            check = 1;
            for n = 1:numel(tmp(:,1))
                ind = find(strcmp(tmp(n,1),ssTotCnt(:,1)));
                tmp(n,2)  = num2cell(cell2mat(tmp(n,2))/cell2mat(ssTotCnt(ind,2)));
            end
            curCnt = {tmp};
        else
            check = 0;
            tmp(1,1) = {''};
            tmp(1,2)  = {0};
            curCnt = {tmp};
        end

    %% scoring

        cnt = [cnt;curCnt];
        curCnt = curCnt{:};
        curCntScore = curCnt;
        curCntScore(:,2) = num2cell(cell2mat(curCntScore(:,2)) * cell2mat(score(i)));
        curCnt(:,2) = num2cell(cell2mat(curCnt(:,2)) * sign(cell2mat(score(i))));


        if check == 1
            curCnt = flipud(sortrows(curCnt,2));
            curCntScore = flipud(sortrows(curCntScore,2));
        end
        tmp = {};
        tmpS ={};    

        for k = 1:size(curCnt,1)
            tmp(k) = strcat(curCnt(k,1),' -> ' ,num2str(curCnt{k,2}));
            tmpS(k) = strcat(curCntScore(k,1),' -> ' ,num2str(curCntScore{k,2}));
        end

        curStr = append(tmp);
        strCnt = [strCnt;curStr];

        curStrScore = append(tmpS);
        strScore = [strScore;curStrScore];

        totCnt = [totCnt;curCnt];
        totScore = [totScore;curCntScore];

    end

    tab = table(paths, cnt,strCnt ,strScore);

    topCount = flipud(sortrows(collapseSum(totCnt),2));
    topScore = flipud(sortrows(collapseSum(totScore),2));


    results_2 = topScore;


    subSystemsPerturbation = sSScoreCompare(model, results1, results_2) ;
else
    subSystemsPerturbation = results1;
end
    























