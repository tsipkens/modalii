
clear;

tic;

t = -100:2:2500; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
% opts.hv = 'constant';
opts.Em = 'default'; %'Krishnan';%'Mie-Krishnan';
prop = Prop({'exper_apb17_Fe','Ar','Fe'}, opts);
% prop = Prop({'exper_apb17_Ag','Ar','Ag'},opts);
% prop = Prop({'exper_ldf','Ar','C'},opts);
prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 298;
prop.Tg = prop.Ti;
prop.sigma = 0;

opts.abs = 'include'; % include absorption
x_fields = {'dp0'};
htmodel = HTModel(prop, x_fields, t, opts);
smodel = SModel(prop, x_fields, t, l);
smodel.htmodel = htmodel;
toc;

tic;
T = htmodel.de_solve((15:15:90)');
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
plot(t, T);
hold on;
plot(t, T_j, 'k--');
hold off;

toc;
