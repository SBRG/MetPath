function [results] = count(what)

% counts the element in a vector

for i=1:numel(what)

    tmpwhat = what(i);
    tmpwhat = tmpwhat{:};
    tosearch=unique(tmpwhat);
    cnt=[];
    for j=1:numel(tosearch)
        cnt(j,1) = sum(strcmp(tosearch(j), tmpwhat));
    end
    if size(tosearch,1)>size(tosearch,2)
        results{i,1} = [tosearch num2cell(cnt)];
    else
        results{i,1} = [tosearch' num2cell(cnt)];
    end
end
        