function newmodel=convertModel(model)

% convert the model in to one compatible with metCHANGE toolbox

model.mets=regexprep(model.mets, '\[', '_');
model.mets=regexprep(model.mets, '\]', '');
model.mets=regexprep(model.mets, '\(', '_');
model.mets=regexprep(model.mets, '\)', '');
model.mets=regexprep(model.mets, '__', '_');
model.mets=regexprep(model.mets, '-', '_');

model.rxns=regexprep(model.rxns, '\[', '_');
model.rxns=regexprep(model.rxns, '\]', '');
model.rxns=regexprep(model.rxns, '\(', '_');
model.rxns=regexprep(model.rxns, '\)', '');
model.rxns=regexprep(model.rxns, '__', '_');
model.rxns=regexprep(model.rxns, '-', '_');
model.rxns=regexprep(model.rxns, '_$', '');

model.grRules = regexprep(model.grRules, '_AT\d{1}' , '');

try
    model.genes = regexprep(model.genes, '_AT\d{1}' , '');
catch
    model.genes = {};
end

newmodel=model;

