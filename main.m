
clear;
clc;
addpath cmap;
tic;

t = -100:2:2500; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
% opts.hv = 'constant';
opts.Em = 'default'; %'Krishnan'; %'Mie-Krishnan';

prop = props.x_apb17_Fe;
prop = props.Ar(prop);

prop = props.Fe(prop, opts);
% prop = props.Ag(prop, opts);
% prop = props.C(prop, opts);

prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 298;
prop.Tg = prop.Ti;
prop.sigma = 0;  % 0.1

% Define models and their parameterizations.
opts.abs = 'include';  % set opts to include absorption
x_fields = {'dp0'};  % set models to take only diameter as inputs
htmodel = HTModel(prop, x_fields, t, opts);
smodel = SModel(prop, x_fields, t, l);
smodel.htmodel = htmodel;
toc;

tic;
T = htmodel.de_solve(prop, (15:15:90)');
toc;

tic;
J = smodel.evaluateF(30);
toc;

tic;
smodel.J = J;
T_j = smodel.evaluateI(30);
toc;

figure(2);
plot(t, squeeze(J));

figure(1);
cmap_sweep(6, flipud(fgreen));  % apply colormap for lines
plot(t, T);
hold on;
plot(t, T_j, 'w:');
hold off;
xlim([min(t), max(t)]);

toc;
