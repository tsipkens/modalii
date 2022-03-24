
% EVALUATE  Primary function to evaluate the heat transfer model at provided x.
%  Acts as a wrapper for DE_SOLVE, updating the values of x.
%  
%  T = HTModel.evaluate(X) evaluates the heat transfer model for a 
%  temperature [K] for the QoI given in X.
%  
%  [T, MPO] = HTModel.evaluate(...) adds an output for the particle mass
%  over time, MPO.
%  
%  [T, MPO, XO] = HTModel.evaluate(...) adds an output for the annealed
%  fraction over time, XO. 
%  
%  AUTHOR: Timothy Sipkens, 2018-11-28

function [T, mpo, Xo] = evaluate(htmodel, x)

%-- Update x values in prop struct ---------------------------------------%
[htmodel, prop] = tools.update_prop(htmodel, x);
%-------------------------------------------------------------------------%

%-- Solve ode using DE_SOLVE function ------------------------------------%
%   Note: Ti is assigned in deSolve directly from the prop structure.
[T, ~, mpo, Xo] = htmodel.de_solve(prop, prop.dp0);
%-------------------------------------------------------------------------%

end

