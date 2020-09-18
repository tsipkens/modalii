
% MAIN_DPG_SQ  Main script demonstrating inference of dpg and sg.
% Timothy Sipkens, 2020-09-17
%=========================================================================%


clear; close all; clc;
addpath('cmap');

t = 0:2:2500; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
opts.Em = 'default';
prop = Prop({'exper_ldf','Ar','C'}, opts);
prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 4150;
prop.Tg = 298;

x_fields = {'dp0','sigma'};
htmodel = HTModel(prop, x_fields, t, opts);
smodel = SModel(prop, x_fields, t, l);
smodel.htmodel = htmodel;



%-- Forward model --------------------------------------------------------%
T = htmodel.de_solve((15:15:90)');

J = smodel.evaluateF([25,log(1.5)]);
J = J ./ max(J(:));

rng(0); % for reproducible results
J_noise = [];
for ii=1:50
    J_noise(:,ii,:) = J + 2e-2 .* sqrt(J) .* randn(size(J)) + 1e-12;
end

% Invert spectroscopy to generate "temperature data".
smodel.J = J_noise;
T_j = smodel.evaluateI([30, log(1.3)]);
%-------------------------------------------------------------------------%



%-- Plot of temperature data ---------------------------------------------%
figure(1);
plot(t, T);
hold on;
plot(t, squeeze(T_j(:,1,:)), 'k--');
hold off;
%-------------------------------------------------------------------------%




%-- Inverse analysis -----------------------------------------------------%
% Setup model
x0 = [22,0.3];
b = @smodel.evaluateI;
model = @smodel.evaluateIF;

% Analysis
tic;
stats = Stats(model,b,opts);
[mle,jcb] = stats.minimize(x0);
disp('MLE = ');
disp(mle);
disp(' ');
toc;

figure(2);
stats.plot_mle(mle);
ylim([0,4500]);


%%
% Evaluate likelihood on a grid for plotting.
% Longer runtimes. 
sigma_vec = linspace(1.0, 2, 78);
dp0_vec = linspace(5, 60, 77);
p1 = zeros(length(dp0_vec), length(sigma_vec));

p_fun = @(x) -(1/2).*norm(stats.min_fun(x)).^2;
disp('Evaluating chi function on grid:');
tools.textbar(0);
for ii=1:length(dp0_vec)
    for jj=1:length(sigma_vec)
        p1(ii, jj) = p_fun([dp0_vec(ii), log(sigma_vec(jj))]);
        tools.textbar((length(sigma_vec) * (ii-1) + jj) ./ ...
            (length(sigma_vec) * length(dp0_vec)));
    end
end
disp('Complete.');
disp(' ');

% Get MCMC samples.
xs = slicesample(mle, 1e3, ...
    'logpdf', @(x) p_fun(x), 'width', 0.1 .* mle);


%%
figure(3);
contourf(sigma_vec, dp0_vec, log(-p1), 40, ...
    'edgecolor', 'none');
colormap(flipud(viridis));
axis square;
set(gca, 'YDir', 'normal');


figure(4);
contourf(sigma_vec, dp0_vec, log(-p1), 40, ...
    'edgecolor', 'none');
colormap(viridis);
axis square;
hold on;
scatter(exp(xs(:,2)), xs(:,1), ...
    'filled', 'w', 'SizeData', 2.5);
alpha(0.5);
hold off;
%-------------------------------------------------------------------------%


%%
%-- Consideration of nuisance parameters ---------------------------------%
y_fields = {'dp0','sigma','Ti','alpha','Tg','Ccp'};
htmodel2 = HTModel(prop, y_fields, t, opts);
smodel2 = SModel(prop, y_fields, t, l);
smodel2.htmodel = htmodel2;
smodel2.J = J_noise;

opts2 = opts;
opts2.minimize = 'vector-xpr'; opts2.prior = 'xpr';

y0 = [22,0.3,4150,0.37,298,1];
sy0 = 0.15.*y0; sy0(1:2) = inf; sy0(3) = 15; sy0(5) = 10;
b2 = @smodel2.evaluateI;
model2 = @smodel2.evaluateIF;
stats2 = Stats(model2,b2,opts2);
stats2.x_pr = y0;
stats2.Lx_ipr = diag(1 ./ sy0);
mle2 = stats2.minimize(y0);

prior_fun = @(x) (x - y0) ./ sy0;
p_fun2 = @(x) -(1/2) .* ( ...
    norm(stats2.min_fun(x)).^2 + ...
    norm(prior_fun(x)).^2);

xs2 = slicesample(mle2, 1e3, ...
    'logpdf', @(x) p_fun2(x), 'width', 0.25 .* mle2);

figure(5);
contourf(sigma_vec, dp0_vec, log(-p1), 40, ...
    'edgecolor', 'none');
colormap(viridis);
axis square;
hold on;
scatter(exp(xs2(:,2)), xs2(:,1), ...
    'filled', 'w', 'SizeData', 2.5);
alpha(0.5);
hold off;
%-------------------------------------------------------------------------%






