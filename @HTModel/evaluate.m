
% EVALUATE  Primary function to evaluate the heat transfer model at provided x.
% Acts as a wrapper for DE_SOLVE, updating the values of x.
% AUTHOR: Timothy Sipkens, 2018-11-28
% 
% INPUTS:
%   x       Vector of numbers with one entry corresponding to each parameter in htmodel.x
%
% OUTPUTS:
%   Tout    Temperature decay for the specified input
%=========================================================================%

function [Tout] = evaluate(htmodel, x)

%-- Update x values in prop struct ---------------------------------------%
[htmodel, prop] = tools.update_prop(htmodel, x);
%-------------------------------------------------------------------------%


%-- Solve ode using DE_SOLVE function ------------------------------------%
%   Note: Ti is assigned in deSolve directly from the prop structure.
Tout = htmodel.de_solve(prop, prop.dp0);
%-------------------------------------------------------------------------%

end

