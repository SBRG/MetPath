function solFinalVals = calculateFluxState(model, mode, allow_loops)

% define flux states
% mode could be: 1-pFBA, 2-low tollerance, 3-minimize fluxes fixing the lb
% of obj, 4- minimize fluxes

if ~exist('allow_loops', 'var')
    allow_loops= 1;
end

if mode == 1 
    if ~exist('model.rules', 'var')
        model.rules = createRules(model);
        display('rules created')
    end
    solFinal = pFBA(model);
    solFinalVals = solFinal.x;
    tol = 10^-6;
    solFinalVals(abs(solFinalVals)<tol) = 0;
end

if mode == 2
    solFinal = optimizeCbModel(model,'max', 10^-6,allow_loops);
    solFinalVals = solFinal.x;
    tol = 10^-6;
    solFinalVals(abs(solFinalVals)<tol) = 0;
end

if mode == 3
    tmp_sol = optimizeCbModel(model);
    model.lb(find(model.c))=tmp_sol.f;
    solFinal = optimizeCbModel(model,'min',10^-6,allow_loops);
    solFinalVals = solFinal.x;
    tol = 10^-6;
    solFinalVals(abs(solFinalVals)<tol) = 0;
end

if mode == 4
    solFinal = optimizeCbModel(model,'min',10^-6,allow_loops);
    solFinalVals = solFinal.x;
    tol = 10^-6;
    solFinalVals(abs(solFinalVals)<tol) = 0;
end
