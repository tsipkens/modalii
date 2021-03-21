
% EVALAUTEI  Bridge function to evaluate spectroscopic inverse model (J > T). 
% AUTHOR: Timothy Sipkens
%=========================================================================%

function [Tout] = evaluateI(smodel, x)

J = smodel.J; % local copy of stored incandescence


%-- Update x values in prop struct ---------------------------------------%
[smodel, prop] = tools.update_prop(smodel, x);
%-------------------------------------------------------------------------%


Tout = smodel.IModel(prop, J); % call IModel function

end

