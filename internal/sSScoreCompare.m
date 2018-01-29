function res = sSScoreCompare(model, ssScore1,ssScore2) 

res = {};
uSS = unique(model.subSystems);
rmv = [];
for i = 1:numel(uSS)
    
    ind1 = find(ismember(ssScore1(:,1), uSS(i)));
    ind2 = find(ismember(ssScore2(:,1), uSS(i)));
    
    if ~isempty(ind1) & ~isempty(ind2)
        res(i,1) = uSS(i);
        res(i,2) = num2cell(log2(cell2mat(ssScore1(ind1,2)) / cell2mat(ssScore2(ind2,2))));
    else
        res(i,1) = uSS(i);
        res(i,2) = {'0'};
        rmv = [rmv;i];
    end
end

try 
    res(rmv,:)=[];
catch
    res(rmv(1:end-1),:)=[];    
end


res = flipud(sortrows(res, 2));
