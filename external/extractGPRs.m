function [parsedGPR,corrRxn] = extractGPRs(model)

warning off all

parsedGPR = [];
corrRxn = [];
cnt = 1;

for i = 1:length(model.rxns)
    if length(model.grRules{i}) > 1
        % Parsing each reactions gpr
        [parsing{1,1},parsing{2,1}] = strtok(model.grRules{i},'or');
        for j = 2:1000
            [parsing{j,1},parsing{j+1,1}] = strtok(parsing{j,1},'or');
            if isempty(parsing{j+1,1})==1
                break
            end
        end
        
        for j = 1:length(parsing)
            for k = 1:1000
                [parsing{j,k},parsing{j,k+1}] = strtok(parsing{j,k},'and');
                if isempty(parsing{j,k+1})==1
                    break
                end
            end
        end
        
        for j = 1:size(parsing,1)
            for k = 1:size(parsing,2)
                parsing{j,k} = strrep(parsing{j,k},'(','');
                parsing{j,k} = strrep(parsing{j,k},')','');
                parsing{j,k} = strrep(parsing{j,k},' ','');
            end
        end
        
        for j = 1:size(parsing,1)-1
            newparsing(j,:) = parsing(j,1:length(parsing(j,:))-1);
        end
        
        parsing = newparsing;
        
     
        for j = 1:size(parsing,1)
            for k = 1:size(parsing,2)
                if length(parsing{j,k}) == 0
                    parsing{j,k} = '';                    
                end
            end
        end
        
              
        num = size(parsing,1);
        for j = 1:num
            sizeP = length(parsing(j,:));
            if sizeP > size(parsedGPR,2)
                for k = 1:size(parsedGPR,1)
                    parsedGPR{k,sizeP} = {''};
                end
            end
            
            for l = 1:sizeP          
            parsedGPR{cnt,l} = parsing(j,l);
            end           
            cnt = cnt+1;
        end
        
        for j = 1:num
            corrRxn = [corrRxn;model.rxns(i)];
        end
        
        clear parsing newparsing
        
    end
    
end

for i = 1:size(parsedGPR,1)
    for j = 1:size(parsedGPR,2)
        if isempty(parsedGPR{i,j}) == 1
            parsedGPR{i,j} = {''};
        end
    end
end

i =1 ;
sizeP = size(parsedGPR,1);
while i < sizeP
    if strcmp(parsedGPR{i,1},{''}) == 1
        parsedGPR = [parsedGPR(1:i-1,:);parsedGPR(i+1:end,:)];
        corrRxn = [corrRxn(1:i-1,:);corrRxn(i+1:end,:)];
        sizeP = sizeP-1;        
        i=i-1;
    end
    i = i+1;
end

for i = 1:size(parsedGPR,1)
    for j= 1:size(parsedGPR,2)
        parsedGPR2(i,j) = cellstr(parsedGPR{i,j});
    end
end

parsedGPR = parsedGPR2;
