function rules = createRules(model)

% create the field rules from grRules

rules = model.grRules;

for i = 1:numel(rules)
    tmp={''};
    tmp =rules(i);
    tmp = strrep(tmp, '(','( ');
    tmp = strrep(tmp, ')',' )');    
    tmp= regexp(tmp, '\s', 'split');
    tmp = tmp{:};
    for k = 1:numel(tmp)
        rep = find(strcmp(tmp(k), model.genes));
        if numel(rep)==1
            tmp(k)=strcat({'x'},{'('},num2str(rep),{')'});
        elseif numel(rep)>1
            display('check')
            i
        end
    end
    tmp = strrep(tmp, 'or', '|');
    tmp = strrep(tmp, 'and', '&');
    rules(i)={''};
     for k = 1:numel(tmp)
        if k==1
            rules(i) = strcat(rules(i),tmp(k));
        else
            rules(i) = strcat(rules(i),{' '},tmp(k));
        end
    end
                
end
