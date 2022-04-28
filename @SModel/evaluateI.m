
% EVALAUTEI  Bridge function to evaluate spectroscopic inverse model (J > T). 
%  
%  T = SModel.evaluateI(X) uses the QoI, X, to evaluate the model.
%  
%  T = SModel.evaluateI(X, J) adds an optional input for incandescence.
%  Otherwise, uses that stored in SModel.
%  
%  AUTHOR: Timothy Sipkens

function [T, C] = evaluateI(smodel, x, J)

if ~exist('J', 'var'); J = []; end
if isempty(J); J = smodel.J; end % use local copy of stored incandescence

%-- Update x values in prop struct ---------------------------------------%
[smodel, prop] = tools.update_prop(smodel, x);
%-------------------------------------------------------------------------%


[T, ~, C] = smodel.IModel(prop, J); % call IModel function

end

