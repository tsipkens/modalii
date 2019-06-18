function [Tout] = evaluateIF(smodel,x)
% EVALAUTEIF Evaluate forward and inverse model (htmodel -> J -> Teff). 
% Used to incorporate polydispersity into forward model. 

htmodel = smodel.htmodel; % embedded heat transfer model
prop = smodel.prop;

if isempty(htmodel) % check if HTModel is missing and produce an error accordingly
    error('Error: Must include HTModel in SModel to use "evaluateIF" function.');
end

if nargin > 1 % update x values
    if length(x)==length(smodel.x)
        for ii=1:length(smodel.x)
            prop.(smodel.x{ii}) = x(ii);
        end
    else
        warning('The number of entries in smodel.x does not match the input dimension.');
    end
end


%-- Consider size distribution -------------------------------------------%
%   Note: In order for the effective temperature to contain polydispersity
%   effects: (i) the HTModel must be evaluated at a range of nanoparticle
%   diameters, (ii) incnadnescence must be generated for the different size
%   classes, (iii) the inandescence must be integrated over the
%   distribution, and (iv) an effectively temperature must be evaluated
%   from that incandescence. 
if prop.sigma > 0.005 % if distribution width is sufficiently large, include polydispersity
    J = smodel.evaluateF; % evaluate forward model for J (includes poly.)
    Tout = smodel.IModel(J); % evaluate inverse model for T
    
else
    Tout = htmodel.evaluate; % for monodisperse case, compute temperature directly
    
end
%-------------------------------------------------------------------------%

end

