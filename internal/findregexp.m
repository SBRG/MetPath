function result=findregexp(str,expr,lind)

% lind = 1 -> result is as logical index 

switch lind
    case 1
        result=~cellfun(@isempty, regexp(str, char(expr)));
    case 0
        result=str(~cellfun(@isempty, regexp(str, char(expr))));
end
end
