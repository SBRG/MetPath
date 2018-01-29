% SOLVE_MILP Pass MILP to solver
%    [x, fval, exitflag] = solve_milp(f, A, b, Aeq, beq, lb, ub, vartype, sense, options)
%    Select SOLVER by uncommenting.
function [x, fval, exitflag] = solve_milp(f, A, b, Aeq, beq, lb, ub, vartype, sense, options)

% SOLVER = 'GLPK';
% SOLVER = 'SCIP';
% SOLVER = 'CPLEX';
SOLVER = 'Gurobi';

if nargin < 10
%     options = optimset('Display', 'final');
    options.Display = 'final';
end
if nargin < 9
    sense = 1;
end
if nargin < 8
    vartype = char('C' * ones(size(f)));
end
if nargin < 6
    lb = [];
    ub = [];
end
if nargin < 4
    Aeq = [];
    beq = [];
end

switch SOLVER
    case 'GLPK'
        param.tmlim = -1;
        param.save = 0;
        switch options.Display
            case 'off'
                param.msglev = 0;
            case 'iter'
                param.msglev = 3;
            case 'final'
                param.msglev = 2;
            case 'notify'
                param.msglev = 1;
        end
        
        ctype = [ char('U' * ones(size(A, 1), 1)); char ('S' * ones(size(Aeq, 1), 1)) ];
        [x, fval, status] = glpk(f, [A; Aeq], [b; beq], lb, ub, ctype, vartype, sense, param);
        
        if status == 180 || status == 151 || status == 171
            exitflag = 1;
        else
            exitflag = -1;
        end
    case 'SCIP'
        param = struct();        
        switch options.Display
            case 'off'
                param.msglev = 0;
            case 'iter'
                param.msglev = 3;
            case 'final'
                param.msglev = 2;
            case 'notify'
                param.msglev = 1;
        end
        
        ctype = [ char('U' * ones(size(A, 1), 1)); char ('S' * ones(size(Aeq, 1), 1)) ];
        if all(vartype == 'C')
            param.tmlim = -1;
            param.save = 0;
            [x, fval, exitflag] = glpk(f, [A; Aeq], [b; beq], lb, ub, ctype, vartype, sense, param);        
        else
            [x, fval, exitflag] = scip(f, [A; Aeq], [b; beq], lb, ub, ctype, vartype, sense, param);        
        end
    case 'CPLEX'
%         cpoptions = cplexoptimset('Display', options.Display);
%         if ~strcmp(options.Display, 'off')
%             cpoptions = cplexoptimset('Diagnostics', 'on');
%         end
        cpoptions = [];
        
        if all(vartype == 'C')
            [x, fval, exitflag] = cplexlp(sense * f, A, b, Aeq, beq, lb, ub, [], cpoptions);
            fval = sense * fval;
        else
            vartype = reshape(vartype, 1, []);
            [x, fval, exitflag] = cplexmilp(sense * f, A, b, Aeq, beq, [], [], [], lb, ub, vartype, [], cpoptions);
            fval = sense * fval;
        end
    case 'Gurobi'
        goptions = struct();
        goptions.FeasibilityTol = 1e-9;
        goptions.IntFeasTol = 1e-9;
        goptions.timeLimit = 200;
%         switch options.Display
%             case 'off'
%                 goptions.Display = 0;
%                 goptions.DisplayInterval = 0;
%             case 'iter'
%                 goptions.Display = 2;
%                 goptions.DisplayInterval = 5;
%             case 'final'
%                 goptions.Display = 2;
%                 goptions.DisplayInterval = 0;
%             case 'notify'
%                 goptions.Display = 1;
%                 goptions.DisplayInterval = 0;
%         end
        
        contypes = [ char('<' * ones(size(A, 1), 1)); 
        char ('=' * ones(size(Aeq, 1), 1)) ];
        
        MILPproblem.lb = lb;
        MILPproblem.ub = ub;
%         MILPproblem.vtype = vartype;
        MILPproblem.vartype = vartype;
%         MILPproblem.obj = f;
        MILPproblem.c = f;
%         MILPproblem.osense = sense;
        MILPproblem.osense = sense;
%         MILPproblem.rhs = [b; beq];
        MILPproblem.b = [b; beq];
%         MILPproblem.b_L = [b; beq];
%         MILPproblem.b_U = [b; beq];
        MILPproblem.A = [A; Aeq];
%         MILPproblem.sense = contypes;
        MILPproblem.csense = contypes;
%         resultgurobi = gurobi(MILPproblem,goptions);
        MILPproblem.x0 = zeros(length(vartype),1);
        resultgurobi = solveCobraMILP(MILPproblem,goptions);
%         [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
%     glnodes, iis, sa, unbdray, farkas] = ...
%     gurobi(f, [A;Aeq], lb, ub, [b; beq], [b; beq], ...
%     'PriLev', 0, 'IntVars', vartype,SC, SI, ...
%     sos1, sos2, logfile, savefile, ...
%     iisRequest, iisFile, saRequest, basis, xIP)
%         resultgurobi
%         stat = resultgurobi.status;
        if resultgurobi.stat==1
%             exitflag = 1;
%             [x, fval, flag] = deal(resultgurobi.x,resultgurobi.objval,resultgurobi.status);
            [x, fval, exitflag] = deal(resultgurobi.full,resultgurobi.obj,resultgurobi.stat);
        else
            x= 0;
            fval = 0;
            exitflag = -1;
        end
        
%         if strcmp(resultgurobi.status,'UNBOUNDED')
%             exitflag = 1;
%         else
%             exitflag = -1;
%         end
end

end
