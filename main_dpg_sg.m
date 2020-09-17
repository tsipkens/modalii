
% MAIN_DPG_SQ  Main script demonstrating inference of dpg and sg.
% Timothy Sipkens, 2020-09-17
%=========================================================================%


clear; close all; clc;

t = -100:2:2500; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
opts.Em = 'default';
prop = Prop({'exper_apb17_Fe','Ar','Fe'}, opts);
prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 298;
prop.Tg = prop.Ti;

opts.abs = 'include'; % include absorption
x_fields = {'dp0','sigma'};
htmodel = HTModel(prop, x_fields, t, opts);
smodel = SModel(prop, x_fields, t, l);
smodel.htmodel = htmodel;



%-- Forward model --------------------------------------------------------%
T = htmodel.de_solve((15:15:90)');

J = smodel.evaluateF([30,log(1.5)]);
J = J ./ max(J(:));

rng(0); % for reproducible results
J_noise = J + 1e-2 .* sqrt(J) .* randn(size(J)) + 1e-12;

% Invert spectroscopy to generate "temperature data".
smodel.J = J_noise;
T_j = smodel.evaluateI([30,log(1.5)]);
%-------------------------------------------------------------------------%



%-- Plot of temperature data ---------------------------------------------%
figure(1);
plot(t, T);
hold on;
plot(t, T_j, 'k--');
hold off;
%-------------------------------------------------------------------------%




%-- Inverse analysis -----------------------------------------------------%
% Setup model
x0 = [35,0.18];
b = @smodel.evaluateI;
model = @smodel.evaluateIF;

% Analysis
tic;
stats = Stats(model,b,opts);
stats.Lb = diag(1./ (100 .* ones(size(T_j))));
[mle,jcb] = stats.minimize(x0);
disp('MLE = ');
disp(mle);
toc;
%-------------------------------------------------------------------------%



