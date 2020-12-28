
% EVALAUTEI  Bridge function to evaluate spectroscopic inverse model (J > T). 
% AUTHOR: Timothy Sipkens
%=========================================================================%

function [Tout] = evaluateI(smodel, x)

J = smodel.J; % local copy of stored incandescence
prop = smodel.prop; % local copy of matl/experimental properties

if nargin > 1 % check if there is a mismatch in size of x
    if length(x)<length(smodel.x)
        error('Error: QoIs parameter size mismatch.');
    else
        if length(x)>length(smodel.x)
            warning('QoIs parameter size mismatch.');
        end
        for ii=1:length(smodel.x)
            prop.(smodel.x{ii}) = x(ii);
        end
    end
end

Tout = smodel.IModel(prop, J); % call IModel function

end

