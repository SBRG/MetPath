function output = escherPaths(met,direction,model,fc, cRes1,modelMets1,cRes2,modelMets2)

% Generate a output containing rxns and the corresponding logFC of rxns scores
% to upload in escher. Since it is considered the logFC the resulting 
% directions of the rxns don't reflect the real directions of the fluxs.

% It is possible use this function just with one condition, in this case it will
% results in rxns scores.

% fc must be a struct object with the following:
%   genes: a vector with gene names
%   expression1: gene expression of the firts condition
%   expression2: gene expression of the second condition

switch direction
    case {'both','b'}
        cRes1.pProdString(ismember(modelMets1.metNames,met));
        rxns1 = intersect(model.rxns, strsplit(char(ans), ';'));
        
        pS1 = cRes1.wProdString(ismember(modelMets1.metNames,met));
        try
            pS1 = str2double(split(pS1, ';'));        
        catch
            pS1 = str2double(strsplit(pS1{1}, ';'));        
        end

        cRes1.pDegString(ismember(modelMets1.metNames,met));
        rxns1 = [rxns1; intersect(model.rxns, strsplit(char(ans), ';'))];
        
        dS1 = cRes1.wDegString(ismember(modelMets1.metNames,met));
        try
            dS1 = str2double(split(dS1, ';'));
        catch
            dS1 = str2double(strsplit(dS1{1}, ';'));
        end

        
        try
            ws1 = [pS1;dS1];
        catch
            ws1 = [pS1';dS1'];
        end

        
        if exist('cRes2', 'var')
            cRes2.pProdString(ismember(modelMets2.metNames,met));
            rxns2 = intersect(model.rxns,strsplit(char(ans), ';'));
            cRes2.pDegString(ismember(modelMets2.metNames,met));
            rxns2 = [rxns2 ; intersect(model.rxns, strsplit(char(ans), ';'))];
            rxns = [rxns2; rxns1];
            
            % adding prodScores and degScores
            pS2 = cRes2.wProdString(ismember(modelMets2.metNames,met));
            
            try
                pS2 = str2double(split(pS2, ';'));
            catch
                pS2 = str2double(strsplit(pS2{1}, ';'));
            end

            dS2 = cRes2.wDegString(ismember(modelMets2.metNames,met));
            
            try
                dS2 = str2double(split(dS2, ';'));
            catch
                dS2 = str2double(strsplit(dS2{1}, ';'));
            end
            
            try
                ws2 = [pS2; dS2];
            catch
                ws2 = [pS2'; dS2'];
            end
            
        else
            rxns = rxns1;
        end
        
      
        
        rxns = unique(rxns);
        genes = {};
        map1 = {};
        map2 = {};
        for i = 1:numel(rxns)
            genes(i,1) = findGenesFromRxns(model, rxns(i));
            map1{i} = find(ismember(rxns1, rxns(i)));
            if exist('cRes2', 'var')
                map2{i} = find(ismember(rxns2, rxns(i)));
            end
        end
        
        if exist('cRes2', 'var')
            for i = 1:numel(genes)
                curLind = ismember(fc.genes, genes{i});
                if sum([~isempty(map1{i}); ~isempty(map2{i})]) == 2
                    rap(i,1) = log((nanmean(fc.expression1(curLind)*(ws1(map1{i}))))./(nanmean(fc.expression2(curLind)*(ws2(map2{i})))));
                elseif ~isempty(map1{i})
                    rap(i,1) = log(nanmean(fc.expression1(curLind)*(ws1(map1{i}))));
                    if rap(i,1) < 0
                        rap(i,1) = rap(i,1)*-1;
                    end                     
                elseif ~isempty(map2{i})
                    rap(i,1) = log(nanmean(fc.expression2(curLind)*(ws2(map2{i}))));
                    if rap(i,1) > 0
                        rap(i,1) = rap(i,1)*-1;
                    end
                end
            end
        else
            for i = 1:numel(genes)
                curLind = ismember(fc.genes, genes{i});
                rap(i,1) = log(nanmean(fc.expression1(curLind)*(ws1(map1{i}))));
            end
        end
        
        
        output = char([{'rxns,val'};strcat(rxns,',', num2str(rap))]);
        

        
    case {'production','p'}
       
        cRes1.pProdString(ismember(modelMets1.metNames,met));
        rxns1 = intersect(model.rxns, strsplit(char(ans), ';'));
        pS1 = cRes1.wProdString(ismember(modelMets1.metNames,met));
        
        try
            pS1 = str2double(split(pS1, ';'));
        catch
            pS1 = str2double(strsplit(pS1{1}, ';'));
        end
 
        ws1 = pS1;
        
        if exist('cRes2', 'var')
            cRes2.pProdString(ismember(modelMets2.metNames,met));
            rxns2 = intersect(model.rxns,strsplit(char(ans), ';'));
            % adding prodScores and degScores
            pS2 = cRes2.wProdString(ismember(modelMets2.metNames,met));
            
            try
                pS2 = str2double(split(pS2, ';'));
            catch
                pS2 = str2double(strsplit(pS2{1}, ';'));
            end
            
            ws2 = pS2;

            rxns = [rxns2; rxns1];
        else
            rxns = rxns1;
        end

        rxns = unique(rxns);
        genes = {};
        map1 = {};
        map2 = {};
        for i = 1:numel(rxns)
            genes(i,1) = findGenesFromRxns(model, rxns(i));
            map1{i} = find(ismember(rxns1, rxns(i)));
            if exist('cRes2', 'var')
                map2{i} = find(ismember(rxns2, rxns(i)));
            end
        end        
        if exist('cRes2', 'var')
            for i = 1:numel(genes)
                curLind = ismember(fc.genes, genes{i});
                if sum([~isempty(map1{i}); ~isempty(map2{i})]) == 2
                    rap(i,1) = log((nanmean(fc.expression1(curLind)*(ws1(map1{i}))))./(nanmean(fc.expression2(curLind)*(ws2(map2{i})))));
                elseif ~isempty(map1{i})
                    rap(i,1) = log(nanmean(fc.expression1(curLind)*(ws1(map1{i}))));
                    if rap(i,1) < 0
                        rap(i,1) = rap(i,1)*-1;
                    end                     
                elseif ~isempty(map2{i})
                    rap(i,1) = log(nanmean(fc.expression2(curLind)*(ws2(map2{i}))));
                    if rap(i,1) > 0
                        rap(i,1) = rap(i,1)*-1;
                    end
                end
            end
        else
            for i = 1:numel(genes)
                curLind = ismember(fc.genes, genes{i});
                rap(i,1) = log(nanmean(fc.expression1(curLind)*(ws1(map1{i}))));
            end
        end
        
        output = char([{'rxns,val'};strcat(rxns,',', num2str(rap))]);
        
    case {'degradation','d'}
        cRes1.pDegString(ismember(modelMets1.metNames,met));
        rxns1 = intersect(model.rxns, strsplit(char(ans), ';'));
        dS1 = cRes1.wDegString(ismember(modelMets1.metNames,met));
        
        try
            dS1 = str2double(split(dS1, ';'));
        catch
            dS1 = str2double(strsplit(dS1{1}, ';'));
        end
                    
        ws1 = dS1;
        
        if exist('cRes2', 'var')
            cRes2.pDegString(ismember(modelMets2.metNames,met));
            rxns2 = intersect(model.rxns,strsplit(char(ans), ';'));
            % adding prodScores and degScores
            dS2 = cRes2.wDegString(ismember(modelMets2.metNames,met));

            try
                dS2 = str2double(split(dS2, ';'));
            catch
                dS2 = str2double(strsplit(dS2{1}, ';'));
            end
            
            ws2 = dS2;

            rxns = [rxns2; rxns1];
        else
            rxns = rxns1;
        end

        rxns = unique(rxns);
        genes = {};
        map1 = {};
        map2 = {};
        for i = 1:numel(rxns)
            genes(i,1) = findGenesFromRxns(model, rxns(i));
            map1{i} = find(ismember(rxns1, rxns(i)));
            if exist('cRes2', 'var')
                map2{i} = find(ismember(rxns2, rxns(i)));
            end
        end        
        if exist('cRes2', 'var')
            for i = 1:numel(genes)
                curLind = ismember(fc.genes, genes{i});
                if sum([~isempty(map1{i}); ~isempty(map2{i})]) == 2
                    rap(i,1) = log((nanmean(fc.expression1(curLind)*(ws1(map1{i}))))./(nanmean(fc.expression2(curLind)*(ws2(map2{i})))));
                elseif ~isempty(map1{i})
                    rap(i,1) = log(nanmean(fc.expression1(curLind)*(ws1(map1{i}))));
                    if rap(i,1) < 0
                        rap(i,1) = rap(i,1)*-1;
                    end                     
                elseif ~isempty(map2{i})
                    rap(i,1) = log(nanmean(fc.expression2(curLind)*(ws2(map2{i}))));
                    if rap(i,1) > 0
                        rap(i,1) = rap(i,1)*-1;
                    end
                end
            end
        else
            for i = 1:numel(genes)
                curLind = ismember(fc.genes, genes{i});
                rap(i,1) = log(nanmean(fc.expression1(curLind)*(ws1(map1{i}))));
            end
        end
        
        output = char([{'rxns,val'};strcat(rxns,',', num2str(rap))]);
        
    otherwise
        display('select direction: both, degradation  or production')
end



printRxnFormula(model, rxns);

end





