
opts.deMethod = 'default';
opts.bFun = 1;
opts.abs = 'include';
opts.pyrometry = 'default';

import('CENIDE.*');
prop = Prop({'CENIDE.experGulder.mat','C.mat','N2.mat'},opts);
prop.F0 = 0.05;
prop.Ti = prop.Tg;

dt = 5;
t = [0:dt:3500]';
l = [500];
x_fields = {};
htmodel = HTModel(prop,x_fields,t,opts);
smodel = SModel(prop,x_fields,t,l,opts);

[Texact,~,m_exact] = htmodel.deSolve(prop.dp0);
Cpre = 500.*m_exact./m_exact(1); % Sublimation effect
J = bsxfun(@times,Cpre,smodel.FModel(Texact,prop.Em));
Jmax = 100;
J = J./max(J).*Jmax;

