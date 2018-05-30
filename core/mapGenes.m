function fMap = mapGenes(model, parsedGPR, corrRxn, genes, exprs1, exprs2)

% map the fold change of genes expression onto active rxns
% the first expression data (exprs1) must be referred to the model
% conditions


FC = exprs1./exprs2;

map.corrRxnUnique = unique(corrRxn);
genesCorr = cell(length(map.corrRxnUnique),1);
map.foldChangeRxn = NaN*ones(length(map.corrRxnUnique),1);
map.foldChangesAll = cell(length(map.corrRxnUnique),1);
genesMap = cell(length(map.corrRxnUnique),1);

for i = 1:length(genesCorr)
    curInds = find(strcmp(map.corrRxnUnique{i},corrRxn));
    curGenes = {};
    
    for j = 1:length(curInds)
        curGenes = union(curGenes,parsedGPR(curInds(j),:));
    end
    
    curBlank = find(cell2mat(cellfun(@isempty,curGenes,'UniformOutput',false)));
    curGenes(curBlank) = [];
    genesCorr{i} = curGenes;
    
    curScores = [];
    curGenesMap = {};
    
    for j = 1:length(curGenes)
        curGene = curGenes{j};
        curInd = find(ismember(genes, curGene));
        if ~isempty(curInd)
            curFold = FC(curInd);
            if curFold == inf
                curFold = exprs1(curInd);
            elseif curFold == nan
                curFold = 0;
            end
            curScores = [curScores;curFold];
            curGenesMap = [curGenesMap;curGene];
        end
    end
    
    if ~isempty(curScores)
        map.foldChangeRxn(i) = nanmean(curScores);
        map.foldChangesAll{i} = curScores;
        map.genesMap{i} = curGenesMap;
    end
    
end

% filter the map pulling out just the reactions with data mapped

rxnIndsData = find(~isnan(map.foldChangeRxn));
fMap.corrRxnUniqueData = map.corrRxnUnique(rxnIndsData);
fMap.foldChangeRxnData = map.foldChangeRxn(rxnIndsData);
fMap.foldChangesAllData = map.foldChangesAll(rxnIndsData);
fMap.genesMapData = map.genesMap(rxnIndsData);

% add subSystems information to the fMap or map

fMap.corrSubsUniqueData = cell(length(fMap.corrRxnUniqueData),1);
for i = 1:length(fMap.corrSubsUniqueData)
        fMap.corrSubsUniqueData{i} = model.subSystems{find(strcmp(fMap.corrRxnUniqueData{i}, model.rxns))};
end