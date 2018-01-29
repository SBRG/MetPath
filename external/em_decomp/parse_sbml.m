% PARSE_SBML Parse SBML file
%    FBAMODEL = PARSE_SBML(SBMLFNAME) parses the SBML file SBMLFNAME.txt
%    and returns the FBA model contained in it as FBAMODEL.
function fbamodel = parse_sbml(sbmlfname)

model = TranslateSBML(sbmlfname, 1);

fbamodel.nrxn = length(model.reaction);
fbamodel.nmetab = length(model.species);
fbamodel.S = sparse(fbamodel.nmetab, fbamodel.nrxn);    % Stoichiometric matrix
fbamodel.rxns = cell(fbamodel.nrxn, 1);
fbamodel.metabs = cell(fbamodel.nmetab, 1);
fbamodel.f = zeros(fbamodel.nrxn, 1);       % biological objective
fbamodel.g = zeros(fbamodel.nrxn, 1);       % synthetic objective
fbamodel.vmin = zeros(fbamodel.nrxn, 1);
fbamodel.vmax = zeros(fbamodel.nrxn, 1);
fbamodel.present = true(fbamodel.nrxn, 1);  % reactions active in model

boundary = zeros(fbamodel.nmetab, 1);

% metabolites & boundary conditions
for imetab = 1:fbamodel.nmetab
    fbamodel.metabs{imetab} = model.species(imetab).id;
    boundary(imetab) = model.species(imetab).boundaryCondition;
end

% parse reactions
for irxn = 1:fbamodel.nrxn
    fbamodel.rxns{irxn} = model.reaction(irxn).id;
    
    for ireactant = 1:length(model.reaction(irxn).reactant)
        fbamodel.S(strcmp(model.reaction(irxn).reactant(ireactant).species, fbamodel.metabs), irxn) ...
            = -model.reaction(irxn).reactant(ireactant).stoichiometry;
    end
    for iproduct = 1:length(model.reaction(irxn).product)
        fbamodel.S(strcmp(model.reaction(irxn).product(iproduct).species, fbamodel.metabs), irxn) ...
            = model.reaction(irxn).product(iproduct).stoichiometry;
    end
    
    for iparam = 1:length(model.reaction(irxn).kineticLaw.parameter)
        param = model.reaction(irxn).kineticLaw.parameter(iparam);
        if strcmp('LOWER_BOUND', param.id)
            fbamodel.vmin(irxn) = param.value;
        elseif strcmp('UPPER_BOUND', param.id)
            fbamodel.vmax(irxn) = param.value;
        elseif strcmp('OBJECTIVE_COEFFICIENT', param.id)
            fbamodel.f(irxn) = param.value;
        end
    end
end

fbamodel.metabs = fbamodel.metabs(~boundary);
fbamodel.S = fbamodel.S(~boundary, :);
fbamodel.nmetab = length(fbamodel.metabs);
