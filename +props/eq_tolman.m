
function [gamma] = eq_tolman(prop,dp,T)

gamma = prop.gammaT(T) ./ (1 + 4 * (prop.delta .* 1e-9) ./ dp);

end

