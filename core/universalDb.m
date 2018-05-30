function [pathways, perturbationScores, universal_ssScore] = universalDb(model,ssScore, data_matched, listCond,parsedGPR,corrRxn, numPerms)

% this function permits to map expression data of the condition of interest
% (expr1) against the second condition in the universal database. It
% returns pathways: predicted set of rxns that compose a pathway; and
% if ssScore = 1 -> universal_ssScore: perturbation of model subSystems

model;
data_matched;
listCond;
parsedGPR;
corrRxn;
numPerms;
ssScore;



genes = data_matched.genes;
exprs1 = data_matched.rank;
exprs2 = data_matched.standardRank;

universal_ut = {};
universal_paths = {};
universal_ssScore_1 = {};

universal_ut_2 = {};
universal_paths_2 = {};
universal_ssScore_2 = {};

path_names = model.mets;
best_score = ones(length(model.metNames),1)*0;
best_path_rxns = cell(length(model.metNames),1);
best_path_lvl = cell(length(model.metNames),1);
best_cond = ones(length(model.metNames),1);
best_cond_str = cell(length(model.metNames),1);

pValLowProd = ones(length(model.metNames),1)*nan;
pValHighProd = ones(length(model.metNames),1)*nan;
pValLowDeg = ones(length(model.metNames),1)*nan;
pValHighDeg = ones(length(model.metNames),1)*nan;
wProdString = cell(length(model.metNames),1);
wDegString = cell(length(model.metNames),1);
scoreDeg = ones(length(model.metNames),1)*nan;
scoreProd = ones(length(model.metNames),1)*nan;
deg_string = cell(length(model.metNames),1);
prod_string = cell(length(model.metNames),1);






for i = 1:numel(listCond)
    i = 1 %THIS APPEARS TO BE A TESTING LINE
    tic;
    [i, numel(listCond)]
    
    try
        load(char(listCond(i)))
    catch
        warning('condition file needs to be in path')
    end
    
    modelMets.metIndsActive = modelMets.metIndsActiveCarbon;
    
    modelList(i,1) = cellstr(char(listCond(i)));
    
    % first condition
    
    %fMap = mapGenes(model, parsedGPR, corrRxn, genes, exprs1, exprs2);
    %tmp_cRes = calcRes(model, modelMets, fMap, calcPaths, numPerms);
    
    
    %REDO BASED ON THE UPDATED FUNCTIONS
    fMapStd = mapGenes(modelStdAdjNoBM, parsedGPR,corrRxn, exprData.genes, ... 
    exprData.anaerobic, exprData.aerobic);
    cutoffDistance = 1;
    cutoffFraction = 0.00;
    pathsStd = metPath(modelStdAdjNoBM, modelMetsAna, metsCurCofInorg, cutoffDistance,cutoffFraction);
    cResStd = calcRes(modelStdAdjNoBM, modelMetsStd, fMapStd, pathsStd, numPerms);
    
    
    % universal ss score
    
    universal_ut{i} = calcAggregateScores(modelMets, tmp_cRes);
    universal_ut{i}(cell2mat(universal_ut{i}(:,4)) == 0,:) = [];
    tmpuniversal_paths = universal_ut{i};
    universal_paths = [universal_paths; tmpuniversal_paths];
    tmp_cRes.pathwaysProd = calcPaths_dist2001.pathwaysProd;
    tmp_cRes.pathwaysDeg = calcPaths_dist2001.pathwaysDeg;
    universal_ssScoreStruct{i} = subSystemsScores(model, tmp_cRes, modelMets);
    tmpuniversal_ssScore = universal_ssScoreStruct{i};
    universal_ssScore_1 = [universal_ssScore_1;tmpuniversal_ssScore];
    
    
    
    % estimating best path
    
    for p = 1:numel(tmp_cRes.pProdString)
        
        cur_path = modelMets.metNames(p);
        ind_ut = ismember(universal_ut{i}(:,1), cur_path);
        cur_score = universal_ut{i}(ind_ut,4);
        res_pos = ismember(path_names, cur_path);
        
        if sum (res_pos) > 0
            prev_best_score = best_score(res_pos);
            cur_score = cell2mat(cur_score);
            if cur_score > prev_best_score
                best_score(res_pos) = cur_score;
                best_path_rxns(res_pos) = strcat(tmp_cRes.pProdString(p),';',tmp_cRes.pDegString(p));
                best_path_lvl(res_pos) = strcat(tmp_cRes.pLevel(p),';',tmp_cRes.dLevel(p));
                
                pValLowProd(res_pos) = tmp_cRes.pValLowProd(p);
                pValHighProd(res_pos) = tmp_cRes.pValHighProd(p);
                pValLowDeg(res_pos) = tmp_cRes.pValLowDeg(p);
                pValHighDeg(res_pos) = tmp_cRes.pValHighDeg(p);
                wProdString(res_pos) = tmp_cRes.wProdString(p);
                wDegString(res_pos) = tmp_cRes.wDegString(p);
                scoreDeg(res_pos) = tmp_cRes.scoresDeg(p);
                scoreProd(res_pos) = tmp_cRes.scoresProd(p);
                prod_string(res_pos) = tmp_cRes.pProdString(p);
                deg_string(res_pos) = tmp_cRes.pDegString(p);
                
                best_cond(res_pos) = i;
                best_cond_str(res_pos) = cellstr(char(listCond(i)));
                 
            elseif cur_score == prev_best_score
                best_cond_str(res_pos) = strcat(best_cond_str(res_pos),';',char(listCond(i)));
            end
        end
    end
    
    
    if ssScore == 1
        
        % second condition
        fMap_2 = mapGenes(model,parsedGPR,corrRxn,genes,exprs2, exprs1);
        tmp_cRes_2 = calcRes(model, modelMets, fMap_2, calcPaths_dist2001);
        
        % universal ss score
        universal_ut_2{i} = calcAggregateScores(modelMets, tmp_cRes_2);
        universal_ut_2{i}(cell2mat(universal_ut_2{i}(:,4))==0,:) = [];
        tmpuniversal_paths_2 = universal_ut_2{i};
        universal_paths_2 = [universal_paths_2; tmpuniversal_paths_2];
        tmp_cRes_2.pathwaysProd = calcPaths_dist2001.pathwaysProd;
        tmp_cRes_2.pathwaysDeg = calcPaths_dist2001.pathwaysDeg;
        universal_ssScoreStruct_2{i} = subSystemsScores(model, tmp_cRes_2, modelMets);
        tmpuniversal_ssScore_2 = universal_ssScoreStruct_2{i};
        universal_ssScore_2 = [universal_ssScore_2;tmpuniversal_ssScore_2];
    end
    
    time = toc;
    estT = time*(numel(listCond)-i)/60;
    display(['estimate time left: ', num2str(estT)])
    
end


pathways = [path_names, num2cell(best_score), best_path_rxns, best_path_lvl, best_cond_str,prod_string,num2cell(pValLowProd),num2cell(pValHighProd),...
    wProdString,num2cell(scoreProd),deg_string,num2cell(pValLowDeg),num2cell(pValHighDeg),wDegString,num2cell(scoreDeg)];
header = {'pathway names', 'best scores', 'rxns', 'levels','models', 'prod rxns','pValLowProd','pValHighProd','prod weightings','prod scores','deg rxns'...
    'pValLowDeg','pValHighDeg', 'deg weightings', 'deg scores'};

rmv = best_score == 0;
pathways(rmv,:) = [];
pathways = flipud(sortrows(pathways, 2));

pathways = [header;pathways];


perturbationScores = collapseMean([universal_paths(:,1),universal_paths(:,4)]);
perturbationScores = flipud(sortrows(perturbationScores,2));

ssScore1 = collapseMean(universal_ssScore_1);


if ssScore == 1
    ssScore2 = collapseMean(universal_ssScore_2);
    universal_ssScore = sSScoreCompare(model,ssScore1,ssScore2);
else
    universal_ssScore = {};
end



