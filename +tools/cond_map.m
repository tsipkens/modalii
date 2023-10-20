
% COND_MAP  A function mapping out the different types of conduction.
%  
%  Run another script to generate a heat transfer model
%  (e.g., main).
%  
%  AUTHOR: Timothy Sipkens, 2020-01-03

function [] = cond_map(htmodel, prop)

d_vec = logspace(log10(1e0), log10(1e7), 300);
Tp = 2000;

% Evaluate three different modes: 
% free molecular, continuum, and transition (Fuchs).
[q0, Kn] = htmodel.q_cond(prop, Tp, d_vec, 'free-molecular');
q1 = htmodel.q_cond(prop, Tp, d_vec, 'continuum');
q2 = htmodel.q_cond(prop, Tp, d_vec, 'fuchs');

q0 = q0 ./ (pi .* d_vec .* 1e-9);
q1 = q1 ./ (pi .* d_vec .* 1e-9);
q2 = q2 ./ (pi .* d_vec .* 1e-9);

% Get plot parameters.
x_vec = Kn;

figure(gcf);
loglog(x_vec, q0);
hold on;
loglog(x_vec, q1);
loglog(x_vec, q2, 'k--');
hold off;

xlim([min(x_vec), max(x_vec)]);
ylim([10 .^ (floor(log10(min(q2)))), ...
      10 .^ (ceil(log10(2 .* max(q1))))]);


end
