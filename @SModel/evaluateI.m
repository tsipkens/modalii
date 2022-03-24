
% EVALAUTEI  Bridge function to evaluate spectroscopic inverse model (J > T). 
%  
%  T = SModel.evaluateI(X) uses the QoI, X, to evaluate the model.
%  
%  AUTHOR: Timothy Sipkens

function [T] = evaluateI(smodel, x)

J = smodel.J; % local copy of stored incandescence


%-- Update x values in prop struct ---------------------------------------%
[smodel, prop] = tools.update_prop(smodel, x);
%-------------------------------------------------------------------------%


T = smodel.IModel(prop, J); % call IModel function

end

