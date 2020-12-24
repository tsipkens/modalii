
% DRUDE Function to evaluate the Drude model of optical properties.
% dp included for compatibility with other Em functions
% Requires these parameters to be defined: prop.omega_p, prop.tau. 
%=========================================================================%

function [Em,n,k] = eq_drude(prop, lambda)

[p,q] = size(lambda);

nu = prop.c ./ (lambda*1e-9);
omega = nu * 2 * pi;

epsilon1 = ones(p,q) - prop.omega_p .^ 2 * prop.tau .^2 ./ ...
    (ones(p, q) + omega .^ 2 * prop.tau .^ 2);
epsilon2 = prop.omega_p .^ 2 .* prop.tau ./ ...
    (omega .* (ones(p, q) + omega .^ 2 * prop.tau .^ 2));

n = sqrt(0.5 * (epsilon1 + sqrt(epsilon1 .^ 2 + epsilon2 .^ 2)));
k = sqrt(0.5 * (-epsilon1 + sqrt(epsilon1 .^ 2 + epsilon2 .^ 2)));

num = 6 * n .* k;
denom = (n .^ 2 - k .^ 2 + 2 * ones(p, q)) .^ 2 + 4 * n .^ 2 .* k .^ 2;
Em = num ./ denom;

end

