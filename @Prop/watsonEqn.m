function [hv] = watsonEqn(prop,T)

hv = prop.hvA().*((1-T./prop.Tcr).^prop.n);

end

