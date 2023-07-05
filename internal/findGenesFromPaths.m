function [Pgenes, Dgenes, PDgenes] = findGenesFromPaths(cRes, model, modelMets,paths)


% find involved genes in paths from genes-reactions matrix

for i = 1:size(paths.pathwaysDeg,1)
    Drxns = intersect(modelMets.rxnsActive,model.rxns(find(paths.pathwaysDeg(i,:))));
    if ~isempty(Drxns)
        genes = findGenesFromRxns(model, Drxns);
        geneslist={};
        for k = 1:numel(genes)
            tmp = {};
            if numel(genes{k})>0
                tmp = genes{k};
                geneslist = [geneslist; tmp];
            end
        end
        
        Dgenes{i} = geneslist;    
    
    else
        Dgenes{i} = {''};
    end
    
    
end



for i = 1:size(paths.pathwaysProd,1)
    Prxns = intersect(modelMets.rxnsActive,model.rxns(find(paths.pathwaysProd(i,:))));
    if ~isempty(Prxns)
        genes = findGenesFromRxns(model, Prxns);
        geneslist={};
        for k = 1:numel(genes)
            tmp = {};
            if numel(genes{k})>0
                tmp = genes{k};
                geneslist = [geneslist; tmp];
            end
        end
        
        Pgenes{i} = geneslist;    
    
    else
        Pgenes{i} = {''};
    end
    
end

for i = 1:size(paths.pathwaysProd,1)
    PDgenes{i} = [Pgenes{i}; Dgenes{i}]';
end
end



