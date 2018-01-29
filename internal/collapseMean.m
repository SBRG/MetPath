function newlist = collapseMean(list)

% collapse a list aggregating the values
% list should contain names in first col and values on second.

newlist={};
tosearch =unique(list(:,1));


for i =1:length(tosearch)
    curName = tosearch(i);
    ind = ismember(list(:,1), curName);
    newlist(i,1) = tosearch(i);
    newlist(i,2) = num2cell(nanmean(cell2mat(list(ind,2))));
end


