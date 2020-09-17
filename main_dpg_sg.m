
% MAIN_DPG_SQ  Main script demonstrating inference of dpg and sg.
% Timothy Sipkens, 2020-09-17
%=========================================================================%


clear;

t = -100:2:2500; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
opts.Em = 'default';
prop = Prop({'exper_apb17_Fe','Ar','Fe'}, opts);
prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 298;
prop.Tg = prop.Ti;
prop.sigma = 0;

opts.abs = 'include'; % include absorption
x_fields = {'dp0','sigma'};
htmodel = HTModel(prop, x_fields, t, opts);
smodel = SModel(prop, x_fields, t, l);
smodel.htmodel = htmodel;


%-- Forward model --------------------------------------------------------%
T = htmodel.de_solve((15:15:90)');

J = smodel.evaluateF([30,0]);

smodel.J = J;
T_j = smodel.evaluateI([30,0]);

figure(2);
plot(t, squeeze(J));

figure(1);
plot(t, T);
hold on;
plot(t, T_j, 'k--');
hold off;

