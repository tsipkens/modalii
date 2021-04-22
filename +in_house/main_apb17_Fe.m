
close all;
clear;
fprintf('/n/n-------------------------------------------------'); % clc;

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

% prop = Prop({['exper_apb17_',signal.matl],...
%     [signal.gas],[signal.matl]},opts);

prop = props.x_apb17_Fe;
prop = eval(['props.', [signal.gas], '(prop, opts)']);
prop = eval(['props.', [signal.matl], '(prop, opts)']);

prop.Ti = 3200; % initial temperature
prop.l = [442, 716]; % wavelengths
prop.alpha = 0.181; % thermal accommodation coefficient


%-- Setup model ----------------------------------------------------------%
x_fields = {'dp0', 'alpha'};
x0 = [35, 0.18];

htmodel = HTModel(prop, x_fields, signal.t, opts);
smodel = SModel(prop, x_fields,...
    signal.t, signal.l, signal, htmodel, opts);

prop.Ti = data.get_peak_temp(signal,smodel); % only used to get Ti to start
htmodel.prop.Ti = prop.Ti;
smodel.prop.Ti = prop.Ti;
smodel.htmodel.prop.Ti = prop.Ti;

b = @smodel.evaluateI;
model = @smodel.evaluateIF;


%-- Analysis -------------------------------------------------------------%
tic;
[b1, Lb, sb] = stats.setup_b(b, 'default');
like = stats.build_like(model, b1, Lb, opts.minimize);
[mle, jcb] = stats.minimize(x0, like, [], opts);
toc;

figure(2);
stats.plot_mle(mle, model, b1, signal.t, sb);

[G_po,R_po,s_po] = stats.cred_linear(jcb);
tools.mle2table(mle, x_fields, 'SPO', s_po);

%-{
%-- Setup model w/ nuisance parameters -----------------------------------%
prop0 = prop;  % save original struct
prop = props.Props(prop);  % convert to Props

theta_fields = {'dp0','alpha','Arho','Brho','Ccp','hvb','Tcr','Tg','Tb','CEmr'};
theta0 = 0.9.*[mle,...
    prop.Arho,prop.Brho,prop.Ccp,prop.hvb,prop.Tcr,...
    prop.Tg,prop.Tb,prop.CEmr];
sx_pr = [inf,inf,abs(theta0(3:end)).*0.05];
sx_pr(end-1) = sx_pr(end-1)*0.1;
sx_pr(end) = sx_pr(end-1)*2;
theta0 = theta0'; sx_pr = sx_pr';
htmodel_t = HTModel(prop,theta_fields,signal.t,signal,opts);
smodel_t = SModel(prop,theta_fields,signal.t,signal.l,signal,htmodel_t,opts);
bModel_t = @smodel_t.evaluateI;
AModel_t = @smodel_t.evaluateIF;
opts.minimize = 'vector-xpr';

tic;
[b1_t, Lb_t, sb_t] = stats.setup_b(bModel_t, 'default');
like_t = stats.build_like(AModel_t, b1_t, Lb_t, opts.minimize);
pr_t = stats.build_prior(theta0, diag(1 ./ sx_pr), opts.minimize);
[mle_t, jcb_t] = stats.minimize(theta0, like_t, pr_t, opts);
toc;

figure(3);
stats.plot_mle(mle_t, AModel_t, b1_t, signal.t, sb_t);

[G_po_t,R_po_t,s_po_t] = stats.cred_linear(jcb_t);
tools.mle2table(mle_t, theta_fields, ...
    'SPO', s_po_t, 'SPR', sx_pr, 'X0', theta0);
%}



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
