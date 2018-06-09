function fluxes = calculateFluxState(model, tolFlux)

% Define flux states
%
% This code uses FBA to estimate the flux state of the system, with the
% additional criteria that the solution will minimize the sum of fluxes
%
% tolFlux is the lowest value for a non-zero flux, i.e. all values below
% will be truncated

sol = optimizeCbModel(model,'min',10^-6);
fluxes = sol.x;
fluxes(abs(fluxes)<tolFlux) = 0;

