function [ind, matched_list, missing, missingInd] = match(list1, list2, verbose)

% useful if you have 2 list that you want to compare but they are not in the same
% order.  ind return the position of elements in list1 in list2.  
% list2(ind) = list1


if ~exist('verbose', 'var')
    verbose = 0;
end

ind = [];

for i = 1:numel(list1)
    x = find(ismember(list2,list1(i)));
    if numel(x)==1
        ind = [ind; x];
    elseif numel(x)>1
        for j = 1:numel(x)
            ind = [ind; x(j)];
        end
    else
        if verbose == 1
            display('one element not found')
            list1(i)
        end
        ind(i,1) = nan;
    end
end

missingInd = find(isnan(ind));

ind = ind(~isnan(ind));


if ~isempty(ind)
    matched_list = list2(ind);
end

if ~isempty(missingInd)
    missing = list1(missingInd);
end













