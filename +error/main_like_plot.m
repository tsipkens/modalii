
clear;
% close all;
clc;
load('+CENIDE\viridis.mat');

Noise.setup_C;

% nn = 500;
% theta = 0.01; % photoelectric and other efficiencies
% sigma = 0.08; % percentage of max, i.e. ~ 0.1 = 10%
% tau = 0.3; % percent variation, i.e. ~ 0.1 = 10%
% [S,S_E,S_std] = Noise.simulate_noise2(J.*theta,theta,sigma,tau,nn); % generate measured signals
% S_var_var = ((nn-1).*((nn-1).*moment(S',4)-(nn-3).*(moment(S',2).^2))./(nn^3))';

% load('+Noise/like_infer_simulate_NelderMead.mat');
%%
% S_vec = S(:);

x0 = [theta,sigma,tau];

mm = 5000;

S_fun = @(x)Noise.S_fun1(S_E,S_std,S_var_var,x,mm);

% [S_poly,S_poly_var] = polyfit(S_E,S_std.^2,2);
% S_poly_L = chol(inv((inv(S_poly_var.R)*inv(S_poly_var.R)')*S_poly_var.normr^2/S_poly_var.df));
%     % estimate of chol. factorization of the covariance matrix of parameters
% S_fun = @(x)Noise.S_fun2(S_E,S_poly,S_poly_L,x,mm);

disp('Evaluating function at grid points...');
ind_ii = theta.*(0:0.04:2);
ind_jj = [0.075,0.08,0.085];%sigma.*[0.95,1,1.05];
ind_kk = tau.*(0.8:0.008:1.2);
% ind_ii = 0.4:0.02:1.6;
% ind_jj = [0.95,1,1.05];
% ind_kk = 0.2:0.05:1.8;
loglike = []; %zeros(1,length(ind));
for ii=1:length(ind_ii);
    disp(['Processing ii = ',num2str(ii),' of ',num2str(length(ind_ii)),'.']);
    for kk=1:length(ind_kk)
        for jj=1:length(ind_jj);
            loglike(ii,jj,kk) = ...
                S_fun([ind_ii(ii),ind_jj(jj),ind_kk(kk)]);
        end
    end
end
disp('Completed.')
disp(' ')

%%
% like = exp(-0.5.*loglike);
% like = exp(-0.5.*loglike./500);

like = loglike;
% ind_ii = ind_ii(1:1:end);
% ind_kk = ind_kk(1:1:end);
% like = like(1:1:end,:,1:1:end);

sum_like = sum(sum(sum(like)));
max_like = max(max(max(like)));
% like = like./max_like;

z_max = 3e4;%max(max(max(like)));
z_min = 0;%min(min(min(like)));%0;
z_nn = 20;
z_ind = z_min:(z_max/(z_nn-1)):z_max;

figure(2);
for jj=1:length(ind_jj)
    subplot(1,length(ind_jj),jj);
    [~,h] = contourf(ind_kk,ind_ii,squeeze(like(:,jj,:)),z_ind);
    colormap(flipud(cm));
%     colormap(cm);
    colorbar;
    set(gca,'CLim',[z_min,z_max]);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
%     set(h,'LineStyle','none');
%     hold on;
%     plot(x0(3),x0(1),'ok');
%     xlim([0.24,0.36]);
%     ylim([0,0.02]);
%     hold off;
end


% print('foo.tiff','-dtiff','-r600');

