
clear;

t = -100:2:2500;
l = [442,716];

opts = [];
% opts.hv = 'constant';
opts.Em = 'default';%'Krishnan';%'Mie-Krishnan';
prop = Prop({'Ar.m','Fe.m','experFe.m'},opts);
% prop = Prop({'Ar.m','C.m','experLDF.m'});
prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 298;
prop.Tg = prop.Ti;
prop.sigma = 0;

opts.abs = 'include';
x_fields = {'dp0'};
htmodel = HTModel(prop,x_fields,t,opts);
smodel = SModel(prop,x_fields,t,l);
smodel.htmodel = htmodel;

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
plot(t,squeeze(J));

figure(1);
plot(t,T);
hold on;
plot(t,T_j,'k--');
hold off;


