function [Tout] = evaluateI(smodel,x)
% EVALAUTEI Evaluate spectroscopic inverse model (J -> T). 

J = smodel.J; % experimental or stored incandescence
prop = smodel.prop; % update x values

if nargin > 1 % update x values
    if length(x)==length(smodel.x)
        for ii=1:length(smodel.x)
            prop.(smodel.x{ii}) = x(ii);
        end
    else
        disp('Warning: QoIs parameter size mismatch.');
    end
end

Tout = smodel.IModel(J);

end

