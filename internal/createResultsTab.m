function finaltab = createResultsTab(modelMets, cRes, filename)


finaltab = [modelMets.metNames, num2cell(cRes.scoresProd),cRes.pSubSyst, cRes.wProdString, cRes.pProdString,cRes.pLevel,num2cell(cRes.pValLowProd), num2cell(cRes.pValHighProd),...
    num2cell(cRes.scoresDeg), cRes.pSubSyst ,cRes.wDegString,cRes.pDegString,cRes.dLevel,num2cell(cRes.pValLowDeg), num2cell(cRes.pValHighDeg)];

header = {'metNames' 'Production_Scores' 'Prodution_subSystems' 'Production_paths_weightings' 'Production_Paths' 'pRxn_distances' 'prod_pVal_Low' 'prod_pVal_High'...
    'Degradation_Scores' 'Degradation_subSystems' 'Degradation_paths_weightings' 'Degradation_Paths' 'dRxn_distances' 'deg_pVal_Low' 'deg_pVal_High'};

finaltab = [header;finaltab];

if exist('filename','var')
    finaltab = table(modelMets.metNames, num2cell(cRes.scoresProd),cRes.pSubSyst, cRes.wProdString, cRes.pProdString,cRes.pLevel,num2cell(cRes.pValLowProd), num2cell(cRes.pValHighProd),...
        num2cell(cRes.scoresDeg), cRes.pSubSyst ,cRes.wDegString,cRes.pDegString,cRes.dLevel,num2cell(cRes.pValLowDeg), num2cell(cRes.pValHighDeg));
    
    finaltab.Properties.VariableNames = {'metNames' 'Production_Scores' 'Prodution_subSystems' 'Production_paths_weightings' 'Production_Paths' 'pRxn_distances' 'prod_pVal_Low' 'prod_pVal_High'...
        'Degradation_Scores' 'Degradation_subSystems' 'Degradation_paths_weightings' 'Degradation_Paths' 'dRxn_distances' 'deg_pVal_Low' 'deg_pVal_High'};

    writetable(finaltab, strcat(char(filename),'.csv'), 'Delimiter', ',','WriteVariableNames',0);

end