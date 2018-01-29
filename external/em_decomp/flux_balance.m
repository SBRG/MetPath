% FLUX_BALANCE Flux-balance analysis of FBA model
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL) performs a basic flux-balance
%    analysis of the FBA model FBAMODEL and returns a biomass-maximizing
%    flux distribution in the vector V.  The maximimum and minimum 
%    synthetic objective possible in a biomass-maximizing flux distribution
%    is given in FMAX and FMIN, respectively.
%
%    [V, FMAX, FMIN] = FLUX_BALANCE(FBAMODEL, QUIET) performs the analysis
%    and supresses screen output if QUIET is set to true.

function [v, fmax, fmin] = flux_balance(fbamodel, quiet)

if nargin < 2
    quiet = false;
end

nrxn   = fbamodel.nrxn;
nmetab = fbamodel.nmetab;

% setup basic flux-balance analysis
c = fbamodel.f;
lb = fbamodel.vmin;
ub = fbamodel.vmax;
   
Aeq = [ fbamodel.S                                                                                                               
        sparse(1:nnz(~fbamodel.present), find(~fbamodel.present), ones(nnz(~fbamodel.present), 1), nnz(~fbamodel.present), nrxn) ];
beq = [ zeros(nmetab, 1);
        zeros(nnz(~fbamodel.present), 1); ];

vartype = char('C' * ones(nrxn, 1));
            
options = optimset('Display', 'off');
            
[x, vbiomass, exitflag] = solve_milp(c, [], [], Aeq, beq, lb, ub, vartype, -1, options);

if ~quiet
    fprintf('Biomass flux:    %f\n', vbiomass);
end

% calculating maximum and minimum synthetic objective possible
% assuming maximum biomass flux
if exitflag > 0
    v = x(1:nrxn);    
    
    c = fbamodel.g;     % synthetic objective
  
    if any(c)
        % maximum possible synthetic flux
        lb_nobiomass = lb;
        ub_nobiomass = ub;
        lb_nobiomass(fbamodel.f > 0) = 0;
        ub_nobiomass(fbamodel.f < 0) = 0;
        [x, max_vsynth] = solve_milp(c, [], [], Aeq, beq, lb_nobiomass, ub_nobiomass, vartype, -1, options);         

        Aeq = [ Aeq;
                fbamodel.f' ];  % condition to ensure maximum biomass in the max/min calculations
        beq = [ beq;
                vbiomass ];

        [x, fmin] = solve_milp(c, [], [], Aeq, beq, lb, ub, vartype, 1, options);
        [x, fmax, exitflag] = solve_milp(c, [], [], Aeq, beq, lb, ub, vartype, -1, options);
        if exitflag > 0
            v = x(1:nrxn);

            if ~quiet
                fprintf('Synthetic flux:  [%f, %f] of %f\n', fmin, fmax, max_vsynth);
            end
        end
    end
else
    fmin = NaN;     % error
    fmax = NaN;
    v = [];
end

