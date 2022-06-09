
% MAIN_2  An updated general function to demonstrate PGM error model.
%  
%  AUTHOR: Timothy Sipkens, 2022-06-08


clear;
close all;
clc;
addpath('cmap')


%-- Setup ----------------------------------------------------------------%
opts.deMethod = 'default';
opts.bFun = 1;
opts.abs = 'include';
opts.pyrometry = 'default';

prop = props.x_cenide_gulder;
prop = props.N2(prop);
prop = props.C(prop);

% prop = Prop({'CENIDE.experGulder.mat','C.mat','N2.mat'},opts);
prop.F0 = 0.05;
prop.Ti = prop.Tg;

dt = 2;
t = [-100:dt:1500]';
l = [500];
x_fields = {};
htmodel = HTModel(prop,x_fields,t,opts);
smodel = SModel(prop,x_fields,t,l,opts);

[Texact, ~, m_exact] = htmodel.de_solve(prop, prop.dp0);
Cpre = 500 .* m_exact ./ m_exact(1); % Sublimation effect
J = bsxfun(@times, Cpre, ...
    smodel.FModel(prop, Texact, prop.Em));
Jmax = 100;
J = J./max(J).*Jmax;
%-------------------------------------------------------------------------%


%%
Jn = pgm.add_noise(J, 1, 0, 0, 1);

figure(1);
plot(t, J);
hold on;
plot(t, Jn);
hold off;

