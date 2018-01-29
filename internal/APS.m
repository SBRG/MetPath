function aggregatePerturbationScores = APS(modelMets, cRes, tab, filename)

%APS: Aggregate Perturbation Score
% create a cell array with production scores, degradation scores and aggregate
% perturbation score for each path
% if tab = 1: the output will be in a table instead of in a cell array
% if filename is present, the reults are written in a .csv file in the working directory 

if ~exist('tab', 'var')
    tab = 0;
end


utilization = zeros(length(cRes.scoresDeg),1);

cRes.scoresDeg(isnan(cRes.scoresDeg)) = 0;
cRes.scoresProd(isnan(cRes.scoresProd)) = 0;


for i=1:numel(cRes.scoresProd)
    utilization(i) = cRes.scoresProd(i) + cRes.scoresDeg(i);
end



if tab ==1
    indup = utilization > 0;
    utilization = num2cell(utilization);
    tmpUtilizationTab = table(modelMets.metNames, num2cell(cRes.scoresProd),num2cell(cRes.scoresDeg),utilization);
    aggregatePerturbationScores = sortrows(tmpUtilizationTab, 4);
    aggregatePerturbationScores.Properties.VariableNames = {'metNames', 'Production_Scores', 'Degradation_Scores', 'Score'};
else

    indup = utilization > 0;
    utilization = num2cell(utilization);
    tmpUtilizationTab = [modelMets.metNames, num2cell(cRes.scoresProd),num2cell(cRes.scoresDeg),utilization];
    aggregatePerturbationScores = sortrows(tmpUtilizationTab, 4);
end

if exist('filename', 'var')
    writetable(aggregatePerturbationScores, strcat(char(filename),'.csv'), 'Delimiter', ',','WriteVariableNames',0);
end
