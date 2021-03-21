
opts.deMethod = 'default';
opts.bFun = 1;
opts.abs = 'include';
opts.pyrometry = 'default';

import('CENIDE.*');
prop = Prop({'inHouse.experFe.mat','Fe.mat','Ar.mat'},opts);
prop.F0 = 0.05;
prop.Ti = prop.Tg;

dt = 0.5;
t = [0:dt:750]';
l = [442];
x_fields = {};
htmodel = HTModel(prop,x_fields,t,opts);
smodel = SModel(prop,x_fields,t,l,opts);

[Texact,~,m_exact] = htmodel.deSolve(prop.dp0);
Cpre = 500.*m_exact./m_exact(1); % Sublimation effect
J = bsxfun(@times,Cpre,smodel.FModel(Texact,prop.Em));
Jmax = 1.7171e+07;%100;
J = J./max(J).*Jmax;

