
% close all;
clear;
clc;

opts.variance = 'variance';
opts.deMethod = 'default';
opts.minimize = 'vector';
opts.estimator = 'likelihood';
opts.display = 'default';
opts.bFun = 1;
opts.Em = 'default';

% Load data ***************************************************************
import('inHouse.*');
load('+inHouse\MoAr_sig1-delay.mat');
% signal.data = signal.data(90:end,:,:);
% signal.t = signal.t(90:end);
prop = Prop({['inHouse.exper',signal.matl,'.m'],...
    ['',signal.gas,'.mat'],['',signal.matl,'.mat']},opts);
prop.l = [442,716];

% Model *******************************************************************
x_fields = {'dp0','sigma'};
x0 = [50,0.2];

htmodel = HTModel(prop,x_fields,signal.t,opts);
smodel = SModel(prop,x_fields,signal.t,signal.l,signal,htmodel);
% prop.Ti = signal.getPeakTemp(smodel); % only used to get Ti to start
prop.Ti = 2510; % temp. set peak temperature
bModel = @smodel.evaluateI;
AModel = @smodel.evaluateIF;

% Analysis ****************************************************************
tic;
stats = Stats(AModel,bModel,prop,opts);
[mle,jcb] = stats.minimize(x0,opts);
disp('MLE = ');
disp(mle);
toc;

% Setup models with nuisance parameters ***********************
theta_fields = {'Arho','Brho','Ccp','hvb','Tcr','Tg','Tb','CEmr'};
theta0 = [prop.Arho,prop.Brho,prop.Ccp,prop.hvb,prop.Tcr,prop.Tg,prop.Tb,prop.CEmr];
st_pr = abs(theta0).*0.05;
st_pr(end-1) = st_pr(end-1)*0.1;
st_pr(end) = st_pr(end-1)*2;
Lt_ipr = sparse(chol(diag(1./(st_pr.^2))));

figure(1);
stats.plotmle(mle,signal.t);
[~,R_po,s_po] = stats.credLinear(jcb);
%{
jcb_theta = stats.jcbEst(mle,theta0);
[G_theta,R_theta,s_theta] = stats.credLinear(jcb_theta);

%{
% [X,Y,x,y] = stats.genGrid(mle,s_po,20);
x = 20:1:75;
y = 0.05:0.01:0.45;
[X,Y] = meshgrid(x,y);

% Plot linear estimate of poterior **************
post_est = mvnpdf([X(:),Y(:)],mle,G_po);
post_est = reshape(post_est,length(y),length(x));
figure;
contourf(x,y,log(post_est),20);
%}
% Plot actual posterior *****************
post = stats.evalPost(X,Y);
figure;
contourf(x,y,post,20);
%}
%{
%{
% MCMC sampling of the posterior ***************
[smp_mcmc,G_mcmc,R_mcmc,s_mcmc] = stats.credMCMC(8000,mle);
hold on;
plot(rnd_smp(:,1),rnd_smp(:,2),'.k');
hold off;
%}

% Bootstrap sampling of the posterior **********
[smp_boot,G_boot,R_boot,s_boot] = stats.credMCMC(8000,mle);
hold on;
plot(smp_boot(:,1),smp_boot(:,2),'.k');
hold off;
%}
