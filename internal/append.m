function string = append(cells)

%beta. works with cells array of strings or double but not with double
%vectors

try
    check = strcmp(class(cells(1,1)), {'double'});
catch
    check = 0;
end

string = {''};

if check == 1
    for k = 1:numel(cells)
        if k==1
            string = strcat(string,num2str(cells(k)));
        else
            string = strcat(string,{';'},num2str(cells(k)));
        end
    end
else
    for k = 1:numel(cells)
        if k==1
            string = strcat(string,cells(k));
        else
            string = strcat(string,{';'},cells(k));
        end
    end
end