
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


%-- Load data -------------------------------%
import('in_house.*');
load('+in_house\MoAr_sig1-delay.mat');
% signal.data = signal.data(90:end,:,:);
% signal.t = signal.t(90:end);
% prop = Prop({['exper_apb17_',signal.matl],...
%     [signal.gas],[signal.matl]},opts);

prop = props.exper_apb17_Mo;
prop = eval(['props.', [signal.gas], '(prop, opts)']);
prop = eval(['props.', [signal.matl], '(prop, opts)']);

prop.l = [442,716];


%-- Model -------------------------------%
x_fields = {'dp0','sigma'};
x0 = [50, 0.2];

htmodel = HTModel(prop, x_fields, signal.t, opts);
smodel = SModel(prop, x_fields,...
    signal.t, signal.l, signal, htmodel, opts);
% prop.Ti = data.get_peak_temp(signal,smodel);  % only used to get Ti to start

prop.Ti = 2510; % temp. set peak temperature
htmodel.prop.Ti = prop.Ti;
smodel.prop.Ti = prop.Ti;
smodel.htmodel.prop.Ti = prop.Ti;

b = @smodel.evaluateI;
A = @smodel.evaluateIF;


%-- Analysis -------------------------------%
tic;
stats = Stats(A,b,prop,opts);
[mle,jcb] = stats.minimize(x0,opts);
toc;


figure(2);
stats.plot_mle(mle, model, b1, signal.t, sb);
[G_po,R_po,s_po] = stats.cred_linear(jcb);

tools.mle2table(mle, x_fields, 'SPO', s_po);

