
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
prop = htmodel.prop;
if nargin > 1
    for ii=1:length(htmodel.x)
        prop.(htmodel.x{ii}) = x(ii);
    end
    if ~(length(x)==length(htmodel.x))
        disp('Warning: QoIs parameter size mismatch in htmodel.evaluate.');
    end
end
%-------------------------------------------------------------------------%


%-- Solve ode using DE_SOLVE function ------------------------------------%
%   Note: Ti is assigned in deSolve directly from the prop structure.
Tout = htmodel.de_solve(htmodel.prop.dp0);
%-------------------------------------------------------------------------%

end

