
% EVALAUTEIF  Evaluate forward and inverse model (htmodel -> J -> Teff).
%  NOTE: Used to incorporate polydispersity into forward model.
%  
%  T = SModel.evaluateIF(X) uses the QoI, X, to evaluate the model.
%  
%  AUTHOR: Timothy Sipkens

function [Tout] = evaluateIF(smodel, x)

htmodel = smodel.htmodel; % embedded heat transfer model

if isempty(htmodel) % check if HTModel is missing and produce an error accordingly
    error('Error: Must include HTModel in SModel to use "evaluateIF" function.');
end


%-- Update x values in prop struct ---------------------------------------%
[smodel, prop] = tools.update_prop(smodel, x);
%-------------------------------------------------------------------------%


%-- POLYDISPERSE ---------------------------------------------------------%
%   Note: In order for the effective temperature to contain polydispersity
%   effects: (i) the HTModel must be evaluated at a range of nanoparticle
%   diameters, (ii) incnadnescence must be generated for the different size
%   classes, (iii) the inandescence must be integrated over the
%   distribution, and (iv) an effectively temperature must be evaluated
%   from that incandescence.
if prop.sigma > 0.005 % if distribution width is sufficiently large, include polydispersity
    J = smodel.evaluateF(x); % evaluate forward model for J (includes poly.)
    Tout = smodel.IModel(prop, J); % evaluate inverse model for T

%-- MONODISPERSE ---------------------------------------------------------%
else
    Tout = htmodel.evaluate(x); % for monodisperse case, compute temperature directly
    
end
%-------------------------------------------------------------------------%

end

