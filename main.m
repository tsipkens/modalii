
% MAIN  A general script to test models. 
%  
%  AUTHOR: Timothy Sipkens

clear;
clc;
close all;

addpath cmap;
tic;

t = -100:2:2500; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
% opts.hv = 'constant';
opts.Em = 'default'; %'Krishnan'; %'Mie-Krishnan';

prop = props.get({'x_apb17', 'Ar', 'Fe'}, [], opts);

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
disp('Completed setup.');
disp(' ');

disp('Getting temperature decays.');
tic;
T = htmodel.de_solve(prop, (15:15:90)');
toc;
disp(' ');

disp('Computing incandescence.');
tic;
J = smodel.evaluateF(30);
toc;
disp(' ');

disp('Inverting sprectroscopic model for temperature.');
tic;
smodel.J = J;
T_j = smodel.evaluateI(30);
toc;
disp(' ');

figure(2);
plot(t, squeeze(J));

figure(1);
cmap_sweep(6, flipud(fmviz));  % apply colormap for lines
plot(t, T);
hold on;
plot(t, T_j, 'w:');
hold off;
xlim([min(t), max(t)]);
