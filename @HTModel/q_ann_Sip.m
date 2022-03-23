
% Q_ANN_SIP Rate of annealing based on the simplified model of Sipkens (2019).
%  
%  AUTHOR: Timothy Sipkens, 2018

function [q,dXdt] = q_ann_Sip(htmodel, prop, T, dp, X)


dp = dp .* 1e-9; % convert to meters so everything is in SI units

E = 4e5; % 5e5 from Newell, J. Appl. Polymer Sci., 1996
k0 = 3e-8 .* 2e13 ./ 5; % 4e12 form Michelsen et al., 2003 prev. 1e5
k = k0 .* exp(-E ./ prop.R ./ T);
    % previous model 1.7e10 and 2e5

DH = 0; % multiply by percentage of defects existing?
    % 1.6e5; % close to sum of the quantities given in Michelsen, J/mol
    % 0.8-1.7e6 J/mol from Michelsen et al., 2003

% dXdt = 3 .* ((1 - X) .^ 2 / 3) ./ dp .* k; % rate of change of fraction of particle
dXdt = 6 .* (1 - X) .^ (2/3) ./ dp .* k;
q = DH .* dXdt .* dp .^ 3 .* pi .* prop.rho(T) ./ 6 ./ prop.M;

end
