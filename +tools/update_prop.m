
% UPDATE_PROP  Update the material property structure using an x.
% AUTHOR: Timothy Sipkens, 2020-12-31
%=========================================================================%

function [obj, prop] = update_prop(obj, x)

prop = obj.prop;
if nargin > 1 % update x values
    if length(x)<length(obj.x)
        error('Error: QoIs parameter size mismatch.');
    else
        if length(x)>length(obj.x)
            warning('QoIs parameter size mismatch.');
        end
        for ii=1:length(obj.x)
            prop.(obj.x{ii}) = x(ii);
        end
    end
end
obj.prop = prop;

end
