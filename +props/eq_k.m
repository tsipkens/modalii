
% EQ_K  Return the thermal conductivity.
%  AUTHOR: Kyle Daun, 2020-12-17
%  MODIFIED: Timothy Sipkens

function [k] = eq_k(T, prop)

coeffs = prop.coeffs;

k = (T<1000) .* (exp(coeffs(3,1).*log(T) + coeffs(3,2)./T + ...
        coeffs(3,3)./T.^2 + coeffs(3,4))) + ...
    (T>=1000) .* (exp(coeffs(4,1).*log(T) + coeffs(4,2)./T + ...
        coeffs(4,3)./T.^2 + coeffs(4,4)));

k = k .* 1e-4;

end

