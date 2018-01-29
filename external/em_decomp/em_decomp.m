% EM_DECOMP Elementary mode decomposition
%    [P, W] = EM_DECOMP(V, FBAMODEL) decomposes the flux distribution V for
%    the FBA model FBAMODEL into a set of elementary modes.  The elementary
%    modes are returned in the n x k matrix P, where n is the length of V
%    and k is the number of elementary modes in the decomposition.  Each
%    elementary mode is stored in a column of P.  The contribution of each
%    elementary mode in P to the flux distribution V is given by the 1 x k
%    weight vector W such that P * W' = V.

function [P, w] = em_decomp(v, fbamodel)

M = 1e2;
EPS = 1e-6;

P = [];
w = [];

% options = optimset('Display', 'off');   % options for Gurobi Solver
options.Display = 'off';

nmetab = length(fbamodel.mets);
nrxn = length(fbamodel.rxns);
S = fbamodel.S;

[y, irxn] = max(abs(v));    % reaction with flux of maximum magnitude

while abs(v(irxn)) > EPS    % stop when we have the zero vector
    % set up the MILP
    % x = [ p q ];
    %   p gives the EM extracted
    %   q is a selection variable for the reactions
    
    % objective function
    c = [ zeros(nrxn, 1);
        ones(nrxn, 1) ];
    
    % equality conditions: Aeq * x = beq
    Aeq = [ S                                                                               sparse(nmetab, nrxn);       % Steady-state
        sparse(1:nnz(v == 0), find(v == 0), ones(nnz(v == 0), 1), nnz(v == 0), nrxn)    sparse(nnz(v == 0), nrxn);  % p(j) = 0 if v(j) = 0
        sparse(1, irxn, 1, 1, nrxn)                                                     sparse(1, nrxn); ];         % chosen p(jk) =/= 0
    beq = [ zeros(nmetab, 1);
        zeros(nnz(v == 0), 1);
        sign(v(irxn)) ];
    
    % inequality conditions: A * x <= b
    A = [ sparse(1:nnz(v), find(v ~= 0), -sign(v(v ~= 0)), nnz(v), nrxn)    sparse(nnz(v), nrxn);       % p conforms to v
        speye(nrxn)                                                       -M * speye(nrxn);           % p(j) is free if q(j) selected
        -speye(nrxn)                                                      -M * speye(nrxn); ];        %
    b = [ zeros(nnz(v), 1);
        zeros(nrxn, 1);
        zeros(nrxn, 1); ];
    
    % bounds
    lb = [ -Inf * ones(nrxn, 1);    % p free
        zeros(nrxn, 1) ];        % q
    ub = [ Inf * ones(nrxn, 1);
        ones(nrxn, 1) ];
    
    % variable types
    vartype = [ repmat('C', nrxn, 1);       % p continuous
        repmat('B', nrxn, 1) ];     % q binary
    
    % solve MILP
    [xopt, fmin, exitflag] = solve_milp(c, A, b, Aeq, beq, lb, ub, vartype, 1, options);
    %     [xopt, fmin, exitflag] = solve_milp(c, A, b, Aeq, beq, lb, ub, vartype, 1);
    
    if exitflag < 0
        %check = 'maybe only one reaction is here'
        p = v;                  % failed to solve - make v the final EM
    else
        p = xopt(1:nrxn);       % extracted EM
    end
    max(abs(S*p));
    p = p / sqrt(sum(p.^2));    % normalise p
    q = xopt((nrxn + 1):end);   % selected reactions
    p((q < 0.5)) = 0;
    wnew = min(v(p ~= 0) ./ p(p ~= 0));     % calculate weighting for maximum contribution of p to v
    max(abs(S*p));
    %    disp_ex(p, fbamodel);
    
    P = [ P, p ];        % store results
    w = [ w wnew ];
    
    v = v - wnew * p;           % remove the contribution of p from v
    max(abs(S*v));
    [y, irxn] = max(abs(v));    % find reaction of those remaining in v with flux of maximum magnitude
end

end
