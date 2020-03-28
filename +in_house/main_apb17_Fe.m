
close all;
clear;
clc;

opts.variance = 'independent';
opts.deMethod = 'default';
opts.minimize = 'vector';
opts.estimator = 'likelihood';
opts.prior = 'independent';
opts.bFun = 1;


%-- Load data ------------------------------------------------------------%
import('in_house.*');
load('+in_house\FeAr_sig17.mat');
% load('+in_house\FeNe_sig72.mat');
% load('+in_house\FeHe_sig43.mat');
prop = Prop({['exper_apb17_',signal.matl],...
    [signal.gas],[signal.matl]},opts);
prop.Ti = 3200; % initial temperature
prop.l = [442,716]; % wavelengths
prop.alpha = 0.181; % thermal accommodation coefficient


%-- Setup model ----------------------------------------------------------%
x_fields = {'dp0','alpha'};
x0 = [35,0.18];

htmodel = HTModel(prop,x_fields,signal.t,opts);
smodel = SModel(prop,x_fields,...
    signal.t,signal.l,signal,htmodel,opts);
prop.Ti = data.get_peak_temp(signal,smodel); % only used to get Ti to start
b = @smodel.evaluateI;
model = @smodel.evaluateIF;


%-- Analysis -------------------------------------------------------------%
tic;
stats = Stats(model,b,opts);
[mle,jcb] = stats.minimize(x0,opts);
disp('MLE = ');
disp(mle);
toc;

figure(2);
stats.plot_mle(mle,signal.t);
[G_po,R_po,s_po] = stats.cred_linear(jcb);


%-- Setup model w/ nuisance parameters -----------------------------------%
theta_fields = {'dp0','alpha','Arho','Brho','Ccp','hvb','Tcr','Tg','Tb','CEmr'};
theta0 = [mle,...
    prop.Arho,prop.Brho,prop.Ccp,prop.hvb,prop.Tcr,...
    prop.Tg,prop.Tb,prop.CEmr];
sx_pr = [inf,inf,abs(theta0(3:end)).*0.05];
sx_pr(end-1) = sx_pr(end-1)*0.1;
sx_pr(end) = sx_pr(end-1)*2;
htmodel_t = HTModel(prop,theta_fields,signal.t,signal,opts);
smodel_t = SModel(prop,theta_fields,signal.t,signal.l,signal,htmodel_t,opts);
bModel_t = @smodel_t.evaluateI;
AModel_t = @smodel_t.evaluateIF;
opts.minimize = 'vector-xpr';
stats_t = Stats(AModel_t,bModel_t,opts,'sx_pr',sx_pr,'x_pr',theta0);
stats_t.get_min_fun;
% [mle_theta,jcb_theta] = stats_t.minimize(theta0);


%-- Uncertainties with nuisance parameters ---------------------------------%
jcb_theta = stats_t.jacob_est(theta0);
[G_theta,R_theta,s_theta] = stats_t.cred_linear(jcb_theta);


%{
[X,Y,x,y] = stats_t.gen_grid(mle,s_po,20);

% Plot linear estimate of poterior
post_est = mvnpdf([X(:),Y(:)],mle,G_po);
post_est = reshape(post_est,length(y),length(x));
figure;
contourf(x,y,log(post_est),20);

% Plot actual posterior
post = stats_t.evalPost(X,Y);
figure;
contourf(x,y,post,20);

%{
% MCMC sampling of the posterior=
[smp_mcmc,G_mcmc,R_mcmc,s_mcmc] = stats_t.credMCMC(8000,mle);
hold on;
plot(rnd_smp(:,1),rnd_smp(:,2),'.k');
hold off;
%}

% Bootstrap sampling of the posterior
[smp_boot,G_boot,R_boot,s_boot] = stats_t.credMCMC(8000,mle);
hold on;
plot(smp_boot(:,1),smp_boot(:,2),'.k');
hold off;
%}
