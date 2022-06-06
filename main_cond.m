
% MAIN_COND  A script demonstrating the different types of conduction.
%  
%  Run another script to generate a heat transfer model
%  s(e.g., in_house.main_apb17_Fe).
%  
%  AUTHOR: Timothy Sipkens, 2020-01-03

d_vec = logspace(log10(1e0), log10(1e7), 300);
Tp = 2000;

[q0, Kn] = htmodel.q_cond(prop, Tp, d_vec, 'free-molecular');
q1 = htmodel.q_cond(prop, Tp, d_vec, 'continuum');
q2 = htmodel.q_cond(prop, Tp, d_vec, 'fuchs');

q0 = q0 ./ (pi .* d_vec .* 1e-9);
q1 = q1 ./ (pi .* d_vec .* 1e-9);
q2 = q2 ./ (pi .* d_vec .* 1e-9);

% Get plot parameters.
addpath cmap; cm = ocean(5);
x_vec = Kn;

figure(1);
loglog(x_vec, q0, 'Color', cm(4, :));
hold on;
loglog(x_vec, q1, 'Color', cm(2, :));
loglog(x_vec, q2, 'k--');
hold off;

xlim([min(x_vec), max(x_vec)]);
ylim([10 .^ (floor(log10(min(q2)))), ...
      10 .^ (ceil(log10(2 .* max(q1))))]);


  
  
%%
prop2 = prop;
p_vec = logspace(4, 11, 100);
q0 = zeros(size(p_vec));
q1 = zeros(size(p_vec));
q2 = zeros(size(p_vec));
Kn = zeros(size(p_vec));

for pp=1:length(p_vec)
    prop2.Pg = p_vec(pp);
    [q0(pp), Kn(pp)] = htmodel.q_cond(prop2, Tp, 100, 'free-molecular');
    q1(pp) = htmodel.q_cond(prop2, Tp, 100, 'continuum');
    q2(pp) = htmodel.q_cond(prop2, Tp, 100, 'fuchs');
end

q0 = q0 ./ (pi .* 100e-9);
q1 = q1 ./ (pi .* 100e-9);
q2 = q2 ./ (pi .* 100e-9);

% Get plot parameters.
addpath cmap; cm = inferno(5);
x_vec = Kn;

figure(2);
loglog(x_vec, q0, 'Color', cm(4, :));
hold on;
loglog(x_vec, q1, 'Color', cm(2, :));
loglog(x_vec, q2, 'k--');
hold off;

xlim([min(x_vec), max(x_vec)]);
ylim([10 .^ (floor(log10(min(q2)))), ...
      10 .^ (ceil(log10(5 .* max(q1))))]);

  
  
  