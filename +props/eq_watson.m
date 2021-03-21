
function [hv] = eq_watson(prop, T)

hv = prop.hvA() .* ((1 - T ./ prop.Tcr) .^ prop.n);

end

