
% GETMU  Return the dynamic viscosity of a gas in units of Ns/m^2. 
%  AUTHOR: Kyle Daun, 2020-12-17
%  MODIFIED: Timothy Sipkens

function [mu] = eq_mu(T, prop)

coeffs = prop.coeffs;

mu = (T<1000) .* (exp(coeffs(1,1).*log(T) + coeffs(1,2)./T + ...
        coeffs(1,3)./T.^2 + coeffs(1,4))) + ...
     (T>=1000) .* (exp(coeffs(2,1).*log(T) + coeffs(2,2)./T + ...
        coeffs(2,3)./T.^2 + coeffs(2,4)));

mu = mu .* 1e-7;

end
