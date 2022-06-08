
% FMODEL  Evaluates forward model (T > J) using given E(m) and temperature decay.
%  
%  JOUT = SModel.FModel(PROP, T, EM) evaluates the spectroscopic model for the given
%  properties (PROP), temperature (T), and absorption function (EM).
%  
%  AUTHOR: Timothy Sipkens, 2017

function [Jo] = FModel(smodel, prop, T, Em, X)

if ~exist('X', 'var'); X = []; end
if isempty(X); X = 0; end

p3 = Em(smodel.l, prop.dp0, X) ./ ...
    (smodel.l .* 1e-9) ./ smodel.data_sc;  % deal with signal scaling
p3 = permute(p3, [1,3,2]);

J = smodel.blackbody(T, smodel.l);  % evaluate Planck's Law

% Multiply Planck's Law by necessary values 
% to return incandescence.
Jo = J .* p3;

end

