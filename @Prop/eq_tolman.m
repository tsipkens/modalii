
function [gamma] = eq_tolman(prop, T, dp)

gamma = prop.gammaT(T)./(1+4*(prop.delta.*1e-9)./dp);

end

