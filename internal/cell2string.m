function string = cell2string(cell, separator)

if ~exist('separator','var')
    separator = ';';
end

string={''};
for k = 1:numel(cell)
        if k==1
            string = strcat(string,cell(k));
        else
            string = strcat(string,separator,cell(k));
        end
    end