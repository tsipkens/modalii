function [S,S_E,S_std,S_ss] = simulate_noise(S_in,theta,gamma,tau,nn)
%SIMULATE_NOISE Simulates Poisson, Gaussian,
%and single-shot noisy signal

mm = length(S_in);
n = randn(1,nn);
S = S_in*ones(1,nn)+...
    tau.*S_in*n+...
    sqrt(theta).*sqrt(S_in*ones(1,nn)+tau.*S_in*n).*randn(mm,nn)+...
    gamma.*randn(mm,nn);

S_std = std(S,[],2);
S_E = mean(S,2);
S_ss = S_in*ones(1,nn)+...
    tau.*S_in*n;

end

