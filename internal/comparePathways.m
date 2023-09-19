function [common, diffPaths] = comparePathways(model1,model2, cRes1,cRes2, modelMets1,modelMets2, paths1, paths2)

aggregatePerturbationScore1 = calcAggregateScores(modelMets1,cRes1);
aggregatePerturbationScore2 = calcAggregateScores(modelMets2,cRes2);

vals1 = cell2mat(aggregatePerturbationScore1(:,4));
vals2 = cell2mat(aggregatePerturbationScore2(:,4));


% comparing paths / different conditions 

[commonPaths, indCommon1, indCommon2] = intersect(modelMets1.metNames, modelMets2.metNames);
[diff1, indDiff1] = setdiff(modelMets1.metNames, modelMets2.metNames);  
[diff2, indDiff2] = setdiff(modelMets2.metNames, modelMets1.metNames); 

%% exclusive paths:

% first condition
for i = 1:numel(indDiff1)
    exclusiveRxnPaths1{i,1} = intersect(modelMets1.rxnsActive, model1.rxns(find(paths1.pathwaysProd(indDiff1(i),:)))); %%
    exclusiveRxnPaths1{i,1} = [exclusiveRxnPaths1{i,1}; intersect(modelMets1.rxnsActive, model1.rxns(find(paths1.pathwaysDeg(indDiff1(i),:))))]; %%
    valInd = ismember(aggregatePerturbationScore1(:,1),diff1(i));
    exclusiveRxnPaths1Scores(i,1) = vals1(valInd);
    exclusiveRxnPaths1str(i,1) = append(exclusiveRxnPaths1{i});
end

% second condition
for i = 1:numel(indDiff2)
    exclusiveRxnPaths2{i,1} = intersect(modelMets2.rxnsActive, model2.rxns(find(paths2.pathwaysProd(indDiff2(i),:)))); %%
    exclusiveRxnPaths2{i,1} = [exclusiveRxnPaths2{i,1}; intersect(modelMets2.rxnsActive, model2.rxns(find(paths2.pathwaysDeg(indDiff2(i),:))))]; %%
    valInd = ismember(aggregatePerturbationScore2(:,1),diff2(i));
    exclusiveRxnPaths2Scores(i,1) = vals2(valInd);
    exclusiveRxnPaths2str(i,1) = append(exclusiveRxnPaths2{i});
end



% sort the exclusive paths by score:   
if exist('exclusiveRxnPaths1Scores')
    [~, indSort1] = sort(exclusiveRxnPaths1Scores,1, 'descend');
    trg1 = 1;
else
    trg1 = 0;
end

if exist('exclusiveRxnPaths2Scores')
    [~, indSort2] = sort(exclusiveRxnPaths2Scores,1, 'descend');
    trg2 = 1;
else
    trg2 =0;
end

% store results

if trg1 == 1 & trg2 ==1
    header = {'paths', 'rxns', 'scores'};
    diffPaths.pathsCondition1 = [diff1, exclusiveRxnPaths1str, num2cell(exclusiveRxnPaths1Scores)];
    diffPaths.pathsCondition1 = diffPaths.pathsCondition1(indSort1,:);
    diffPaths.pathsCondition1 = [header; diffPaths.pathsCondition1];
    diffPaths.pathsCondition2 = [diff2, exclusiveRxnPaths2str, num2cell(exclusiveRxnPaths2Scores)];
    diffPaths.pathsCondition2 = diffPaths.pathsCondition2(indSort2,:);
    diffPaths.pathsCondition2 = [header; diffPaths.pathsCondition2];
else
    diffPaths = {};
end

%% common paths:
% extract rxns from common paths:

% first condition
for i = 1:numel(indCommon1)
    commonPathRxns1{i,1} = model1.rxns(find(paths1.pathwaysProd(indCommon1(i),:)));
    commonPathRxns1{i,1} = [commonPathRxns1{i,1}; model1.rxns(find(paths1.pathwaysDeg(indCommon1(i),:)))];

    % filtering out the inactive rxns
    commonPathRxns1{i,1} = intersect(commonPathRxns1{i,1}, modelMets1.rxnsActive);
    div1(i) = numel(commonPathRxns1{i,1});
end

% second condition
for i = 1:numel(indCommon2)
    commonPathRxns2{i,1} = model2.rxns(find(paths2.pathwaysProd(indCommon2(i),:)));
    commonPathRxns2{i,1} = [commonPathRxns2{i,1}; model2.rxns(find(paths2.pathwaysDeg(indCommon2(i),:)))];

    % filtering out the inactive rxns
    commonPathRxns2{i,1} = intersect(commonPathRxns2{i,1}, modelMets2.rxnsActive);
    div2(i) = numel(commonPathRxns2{i,1});
end


% overlap Score
for i = 1:numel(commonPathRxns2)
    overlapScore(i) = sum(ismember(commonPathRxns2{i}, commonPathRxns1{i}))./((numel(commonPathRxns2{i}) + numel(commonPathRxns1{i}))/2); %%
end


% saving the exclusive rxns in common paths in first and second condition
% and their scores

for i = 1:numel(commonPaths)

    % first
    exclusiveRxns1{i,1} = setdiff(commonPathRxns1{i}, commonPathRxns2{i}); %%
    valInd = ismember(aggregatePerturbationScore1(:,1),commonPaths(i));
    exclusiveRxns1Score(i,1) = vals1(valInd);
    exclusiveRxns1str(i,1) = append(exclusiveRxns1{i});

    % second
    exclusiveRxns2{i,1} = setdiff(commonPathRxns2{i}, commonPathRxns1{i}); %%
    valInd = ismember(aggregatePerturbationScore2(:,1),commonPaths(i));
    exclusiveRxns2Score(i,1) = vals2(valInd);
    exclusiveRxns2str(i,1) = append(exclusiveRxns2{i});
end



% shared rxns in common paths

for i = 1:numel(commonPathRxns1)
    [commonRxns{i,1}, commonRxns1Ind{i,1}, commonRxns2Ind{i,1}] = intersect(commonPathRxns1{i}, commonPathRxns2{i});
    valInd = ismember(aggregatePerturbationScore1(:,1),commonPaths(i));
    commonRxnsScr1(i,1) = vals1(valInd);
    valInd = ismember(aggregatePerturbationScore2(:,1),commonPaths(i));
    commonRxnsScr2(i,1) = vals2(valInd);
end


% scoring the fc

for i =1:numel(commonPathRxns1)
    commonRxnsScrFc(i,1) = abs(log(commonRxnsScr2(i,1)./commonRxnsScr1(i,1))); 
    if ~isnan(commonRxnsScrFc(i,1))
        exclusiveRxns1Scorestr(i,1) = num2cell(commonRxnsScr1(i,1));
        exclusiveRxns2Scorestr(i,1) = num2cell(commonRxnsScr2(i,1));
        commonRxnsSortedstr(i,1) = append(commonRxns{i});
        commonFlxsFcSortedstr(i,1) = append(num2str(commonRxnsScrFc(i,1)));
    else
        exclusiveRxns1Scorestr(i,1) = {['0']};
        exclusiveRxns2Scorestr(i,1) = {['0']};
        commonRxnsSortedstr(i,1) = append('');
        commonFlxsFcSortedstr(i,1) = append('');
    end
end


% index for sorting by overlapScore
[overlapScoreSorted, indOLSS] = sort(overlapScore,2,'descend');

[exclusiveRxns1ScoreSorted, indERF1SS] = sort(exclusiveRxns1Score, 1, 'descend');

[exclusiveRxns2ScoreSorted,indERF2SS] = sort(exclusiveRxns2Score, 1, 'descend');



% num2cell
overlapScoreSorted = num2cell(overlapScoreSorted)';
exclusiveRxns1ScoreSorted = num2cell(exclusiveRxns1ScoreSorted)';
exclusiveRxns2ScoreSorted = num2cell(exclusiveRxns2ScoreSorted)'; 



%% results: 

% not sorted 
header1 = {'metName', 'overlap', 'common_rxns', 'fluxes_fc_in_common_rxns', 'rxns_only_in_first','perturbationScore1', 'rxns_only_in_second', 'perturbationScore2'};
common.results_not_sorted = [commonPaths, num2cell(overlapScore)', commonRxnsSortedstr, commonFlxsFcSortedstr,exclusiveRxns1str,exclusiveRxns1Scorestr,exclusiveRxns2str,exclusiveRxns2Scorestr];
common.results_not_sorted = [header1 ; common.results_not_sorted];

%  sorted by overlap score 
[~, sortInd] = sort(overlapScore,2,'descend');
overlapScore = num2cell(overlapScore)';

common.results_overlasp_score = flipud([commonPaths(sortInd), overlapScore(sortInd), commonRxnsSortedstr(sortInd), commonFlxsFcSortedstr(sortInd), exclusiveRxns1str(sortInd),exclusiveRxns1Scorestr(sortInd),exclusiveRxns2str(sortInd),exclusiveRxns2Scorestr(sortInd)]);
common.results_overlasp_score = [header1 ; common.results_overlasp_score];


%  sorted by perturbation score for first condition
header2 = {'metName', 'rxns_only_in_first', 'perturbationScore1','rxns_only_in_second', 'perturbationScore2','overlap', 'common_rxns', 'fluxes_fc_in_common_rxns'};
common.results_sorted_by_score1 = [commonPaths(indERF1SS), exclusiveRxns1str(indERF1SS),exclusiveRxns1Scorestr(indERF1SS),exclusiveRxns2str(indERF1SS),exclusiveRxns2Scorestr(indERF1SS),overlapScore(indERF1SS),commonRxnsSortedstr(indERF1SS), commonFlxsFcSortedstr(indERF1SS)];
common.results_sorted_by_score1 = [header2; common.results_sorted_by_score1];

%  sorted by perturbation score for second condition

header3 = {'metName', 'rxns_only_in_second', 'perturbationScore2','rxns_only_in_first', 'perturbationScore1','overlap', 'common_rxns', 'fluxes_fc_in_common_rxns'};
common.results_sorted_by_score2 = [commonPaths(indERF2SS), exclusiveRxns2str(indERF2SS),exclusiveRxns2Scorestr(indERF2SS),exclusiveRxns1str(indERF2SS),exclusiveRxns1Scorestr(indERF2SS),overlapScore(indERF2SS),commonRxnsSortedstr(indERF2SS), commonFlxsFcSortedstr(indERF2SS)];
common.results_sorted_by_score2 = [header3; common.results_sorted_by_score2];
