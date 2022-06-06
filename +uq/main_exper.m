
clear;
% close all;
clc;
load('+CENIDE\viridis.mat');

% Using simulated C-N2 data
% load('+Noise\data_simulate.mat');
% [W,S_E,S_var,nn,mm,S_var_var] = Noise.calcWeight2(S);
% S_var_var = max(S_var_var).*ones(size(S_var_var));
% W = spdiags(1./S_var_var,0,mm,mm);
% x0 = [tau^2,theta,gamma^2]';

% CENIDE, C, Gulder burner data
% load('+Noise\data_CENIDE_C.mat');
% S = J4;
% [W,S_E,S_var,nn,mm,S_var_var] = Noise.calcWeight2(S);
% S_var_var = max(S_var_var).*ones(size(S_var_var));
% W = spdiags(1./S_var_var,0,mm,mm);

% NRC, C, 2015
% load('+Noise\data_NRC_C.mat');
% S = J2;
% [W,S_E,S_var,nn,mm,S_var_var] = Noise.calcWeight2(S);
% S_var_var = max(S_var_var).*ones(size(S_var_var));
% W = spdiags(1./S_var_var,0,mm,mm);

% CENIDE, Si, 2014
% load('+Noise\data_CENIDE14_Si.mat');
% S = J1; 
% [W,S_E,S_var,nn,mm,S_var_var] = Noise.calcWeight(S);
% S_var_var = max(S_var_var).*ones(size(S_var_var));
% W = spdiags(1./S_var_var,0,mm,mm);

% CENIDE, Ge, plasma synthesized, PMTs same as laminar diffusion flame
% load('..\Data\Ge (CENIDE)\Ge_2mJ_singleshot_PMTs.mat');
% S = ch1';
% [W,S_E,S_var,nn,mm,S_var_var] = Noise.calcWeight2(S);
% S_var_var = max(S_var_var).*ones(size(S_var_var));
% W = spdiags(1./S_var_var,0,mm,mm);

% In-house, Fe-Ar, aerosolizer data
% load('+inHouse\FeHe_sig43.mat');
% load('+inHouse\FeAr_sig17.mat');
load('+inHouse\FeNe_sig72.mat');
S = squeeze(signal.J_raw(2:900,:,1));
[W,S_E,S_var,nn,mm,S_var_var] = Noise.calcWeight2(S);
S_var_var = max(S_var_var).*ones(size(S_var_var));
W = spdiags(1./S_var_var,0,mm,mm);

% DLR, turbulent soot
% load('+Noise\data_DLR_C.mat');
% S_raw = data(20:160,25:35,:); % used (20:160,25;35) on 6/3/2017
% S_raw = S_raw./(max(max(mean(S_raw,3)))).*100;
% S_sz = size(S_raw);
% S = reshape(S_raw,S_sz(1).*S_sz(2),S_sz(3));
% S = S./max(mean(S,2)).*100;
% S = S(mean(S,2)<35,:);
% [W,S_E,S_var,nn,mm] = Noise.calcWeight(S);

A = [S_E.^2,S_E,ones(length(S_E),1)];
B = (A')*W;
% x_mle = fmincon(@(x) (A*x-S_var)'*W*(A*x-S_var),[0.02;12;0.1],...
%     [],[],[],[],[0;0;0],[inf;inf;inf]);
x_mle = fmincon(@(x) (A*x-S_var)'*W*(A*x-S_var),[0.0002;0.07;0.001],...
    [],[],[],[],[0;0;0],[inf;inf;inf]);
% x_mle = inv(B*A)*B*S_var;
tau = sqrt(x_mle(1));
theta = x_mle(2);
gamma = sqrt(x_mle(3));

% A = [S_E,ones(length(S_E),1)];
% B = (A')*W;
% x_mle = inv(B*A)*B*S_var;
% tau = 0;
% theta = x_mle(1);
% gamma = sqrt(x_mle(2));

% Plot E vs var **********************************************************
disp('Generating figure...');
figure(3);
x_max = 100;
plot(S_E,S_var,'.');
hold on;
S_quad = @(E) (tau^2).*(E.^2)+theta.*E+gamma^2;
fplot(S_quad,[0,x_max]);
fplot(@(E) theta.*E+gamma^2,[0,x_max],'-k');
xlim([0,x_max]);
hold off;
disp('Completed');
disp(' ');

G_post = inv(B*A); % posterior covariance matrix
R_post = corrcov(G_post);
s_tau2 = sqrt(G_post(1,1));
s_theta = sqrt(G_post(2,2));
s_gamma2 = sqrt(G_post(3,3));

mrnd_post = mvnrnd(x_mle,G_post,10000); % samples for uncertainties in tau and gamma
s_tau = std(sqrt(mrnd_post(:,1)));
s_gamma = std(sqrt(mrnd_post(:,3)));

% output to the screen
disp('Maximum likelihood estimates (MLEs): ');
disp(['tau^2 = ',num2str(tau^2),' ',char(177),num2str(s_tau2)]);
disp(['tau = ',num2str(tau),' ',char(177),num2str(s_tau)]);
disp(['theta = ',num2str(theta),' ',char(177),num2str(s_theta)]);
disp(['gamma^2 = ',num2str(gamma^2),' ',char(177),num2str(s_gamma2)]);
disp(['gamma = ',num2str(gamma),' ',char(177),num2str(s_gamma)]);
disp(['sigma_rel = ',num2str(gamma/theta)]);
disp(' ');

S_star = (-theta-sqrt(theta^2-4*gamma^2*(tau^2-1)))/(2*(tau^2-1));
disp('S_star(SNR=1) = ');
disp(['   ',num2str(S_star)]);
disp(' ');

%{
%% Generate contour plots *************************

S_x1 = @(x) x(1).*(S_E.^2)+x(2).*S_E+x(3);
S_x2 = @(x) (x(1)^2).*(S_E.^2)+x(2).*S_E+(x(3)^2);
L = chol(W); % Cholesky factorization of inverse covariance matrix

loglike1_fun = @(x) mvnpdf([x(1),x(2),x(3)],x_mle',G_post); %0.5.*norm(((S_x1(x)-S_var)')*L).^2;
loglike2_fun = @(x) mvnpdf([x(1).^2,x(2),x(3).^2],x_mle',G_post); %0.5.*norm(((S_x2(x)-S_var)')*L).^2;
% ind_diff = 3.*sqrt(diag(G_post))./x0;
% ind_tau = sqrt(x_mle(1).*((1-ind_diff(1)):ind_diff/40:(1+ind_diff(1))));
% ind_theta = x_mle(2).*((1-ind_diff(2)):ind_diff/50:(1+ind_diff(2)));
% ind_gamma = sqrt(x_mle(3).*((1-ind_diff(3)):ind_diff/90:(1+ind_diff(3))));
ind_tau = 0.10:0.00025:0.14;
ind_theta = 0.6:0.0025:1.262;
ind_gamma = 4.1:0.005:4.9;

figure(2);
F = [];
for ii=1:length(ind_tau);
    disp(['ii = ',num2str(ii),' of ',num2str(length(ind_tau))]);
    for jj=1:length(ind_theta);
        for kk=1:length(ind_gamma);
            ind_x = [ind_tau(ii).^2,ind_theta(jj),ind_gamma(kk).^2];
%             F(ii,jj,kk) = -0.5.*((A*ind_x-S_var)')*W*(A*ind_x-S_var);
            F(ii,jj,kk) = mvnpdf(ind_x,x_mle',G_post);
        end
    end
end

like = log(squeeze(sum(F,2)));

% mrnd(:,1) = sqrt(mrnd(:,1));
% mrnd(:,3) = sqrt(mrnd(:,3));

sum_like = sum(sum(sum(like)));
max_like = max(max(max(like)));
% like = like./max_like;
min_like = min(min(min(like)));
% like = like./max_like;

z_max = 0;%max_like;
z_min = -90;%min_like;
z_nn = 8;
z_ind = z_min:((z_max-z_min)/(z_nn)):z_max;

contourf(ind_gamma,ind_tau,like,9);
colormap(cm);
% colorbar;
% set(gca,'CLim',[z_min,z_max]);
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     set(h,'LineStyle','none');

%%
mrnd = mvnrnd(x_mle,G_post,400);
plot(sqrt(mrnd(:,3)),sqrt(mrnd(:,1)),'.');

hold on;
plot(sqrt(x0(3)),sqrt(x0(1)),'ok');
xlim([4.1,4.9]);
ylim([0.1,0.14]);
hold off;

hold on;
plot(sqrt(x_mle(3)),sqrt(x_mle(1)),'sk');
xlim([4.1,4.9]);
ylim([0.1,0.14]);
hold off;

%}
