
% ADD_NOISE  Simulates Poisson, Gaussian, and multiplicative noise.

function [s, se, s_std, s_ss] = add_noise(s_in, theta, gamma, tau, nn)

mm = length(s_in);
n = randn(1,nn);
s = s_in * ones(1,nn) + ...
    tau .* s_in * n + ...
    sqrt(theta) .* sqrt(s_in * ones(1,nn) + tau .* s_in * n) .* randn(mm,nn)+...
    gamma .* randn(mm,nn);

s_std = std(s, [], 2);
se = mean(s, 2);
s_ss = s_in * ones(1, nn) + ...
    tau .* s_in * n;

end

