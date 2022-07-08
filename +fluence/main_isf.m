
clear;
clc;
close all;

addpath cmap;
tic;

t = -20:0.1:100; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
% opts.hv = 'constant';
opts.Em = 'default'; %'Krishnan'; %'Mie-Krishnan';

prop0 = props.x_ldf;
prop0 = props.Ar(prop0);

prop0 = props.C(prop0, opts);

prop0.Ti = 298;  % for cold soot
prop0.Tg = prop0.Ti;

prop0.sigma = 0;  % 0.1

prop = prop0;  % copy to simple model, prior to changes

%-- Change true model parameters -----------------------------------------%
% Case 1: Simple doubling of E(m).
% prop0.Em = @(l,dp,X) 0.8 .* ones(1,length(l));

% Case 2: Annealing.
%{
opts.ann = 'Sipkens';
bet0 = 28.72;  % soot
zet0 = 0.83;
bet1 = 0.626;  % graphite
zet1 = 1.186;
fun_Em = @(l,dp,bet,zet) (l .* 1e-6) .^ (1-zet) .* bet ./ (6*pi);
prop0.Em = @(l,dp,X) (1-X) .* fun_Em(l,dp,bet0,zet0) + X .* fun_Em(l,dp,bet1,zet1);
prop0.Emr = @(l1,l2,dp) prop0.CEmr.*prop0.Em(l1,dp,0)./prop0.Em(l2,dp,0);
%}

% Case 3.
% See loop below.

% Case 4.
%{
% Keep this commented for Tg = 298 K (cold soot) case.
prop0.Ti = 1730;  % fot hot soot
prop0.Tg = prop0.Ti;
prop = prop0;
% + See loop below.
%}

% Case 5.
% See below loop. 
%-------------------------------------------------------------------------%

% Define models and their parameterizations.
opts.abs = 'include';  % set opts to include absorption
x_fields = {'dp0', 'F0', 'CEmr'};  % set models to take only diameter as inputs

htmodel = HTModel(prop0, x_fields, t, opts);

% True SModel, used for generating signals.
smodel = SModel(prop0, x_fields, t, l);
smodel.htmodel = htmodel;

% Simple SModel.
smodels = SModel(prop, x_fields, t, l);
smodels.htmodel = htmodel;

disp('Completed setup.');
disp(' ');

[~, prompt] = min(abs(t - 20));

nf = 200;
% nf = 20;
F0_vec = linspace(0.0005, 0.35, nf);
dp = 40;

T = [];  J1 = [];  J2 = []; Cinf = []; Ct = [];
disp('Computing temperature decays:');
tools.textbar([0, nf]);
for ii=1:length(F0_vec)
    [T(:,ii), m(:,ii), X(:,ii)] = htmodel.evaluate([dp, F0_vec(ii), 1]);
    
    Jt = smodel.evaluateF([dp, F0_vec(ii), 1]);
    Jt = Jt .* m(:,ii) ./ m(1,ii); % scale by particle mass loss

    % For Case 3: Time response.
    %{
    tg = 2;  % [2,5]
    Jt = sum(permute(Jt, [2,1,3]) .* normpdf(t - t' + tg, 0, tg), 2) ./ ...
        sum(normpdf(t - t', 0, tg), 2);
    %}
    
    % For Case 4: Background luminance.
    %{
    Jt = Jt + smodel.evaluateF([dp, 0, 1]);
    %}
    
    J1(:,ii) = Jt(:,1,1);  % choose first wavlength
    J2(:,ii) = Jt(:,1,2);  % choose second wavlength
    
    [~, ~, Ct(:,ii)] = smodel.IModel(prop0, Jt);  % inferred ISF (not necessarily true)
    
    % Js = smodels.evaluateF([dp, F0_vec(ii), 1]);
    [~, ~, Cinf(:,ii)] = smodels.IModel(prop, Jt);  % inferred ISF
    
    tools.textbar([ii, nf]);
end
disp(' ');
Jn1 = J1;  Jn2 = J2;  % copy over to temporary array, in case modified below


%%
%-{
% For Case 5.
% LII workshop: (0.2, 0.005)
sF0_1 = 0.;
sF0_2 = 0.05;
disp('Convolving fluence responses (simulating uneven fluence):');
tools.textbar([0, nf]);
for ii=1:length(F0_vec)
    sF0 = sF0_1 * F0_vec(ii) + sF0_2;
    if or(or(F0_vec(ii) < 0.01, ...
            (F0_vec(ii) - 2 * sF0) < 0), ...
            (F0_vec(ii) + 2 * sF0) > F0_vec(end))
        Cinf(:,ii) = NaN;

        tools.textbar([ii, nf]);
        continue;
    end

    pf = normpdf(F0_vec, F0_vec(ii), sF0);
    pf = pf ./ sum(pf);  % scale PDF, so that ISF for signals ideally equals one
    
    J1(:,ii) = sum(Jn1 .* pf, 2);
    J2(:,ii) = sum(Jn2 .* pf, 2);
    Jt = cat(3, J1(:,ii), J2(:,ii));
    
    % Js = smodels.evaluateF([dp, F0_vec(ii), 1]);
    [~, ~, Cinf(:,ii)] = smodels.IModel(prop, Jt);  % inferred ISF

    
    tools.textbar([ii, nf]);
end
srel = sF0 ./ F0_vec;
%}



cm = flipud(internet);

figure(1);
cmap_sweep(nf, cm);
plot(t, T);

Tlow = 10000 * 6 * pi * prop0.Eml(dp) / ...
    (prop0.l_laser * 1e-9 * prop0.rho(1e3) * prop0.cp(1e3)) ...
    .* F0_vec + prop0.Tg;

figure(2);
plot(F0_vec, max(T));
hold on;
plot(F0_vec, Tlow)
hold off;

% Peak.
figure(3);
plot(F0_vec, max(J1));
hold on;
plot(F0_vec, max(Jn1 .* m ./ m(1,:)));
hold off;
title('Peak incandescence v. fluence');

% Constant time from peak laser energy.
np = 400;
figure(4);
plot(F0_vec, J1(np,:));
hold on;
plot(F0_vec, (Jn1(np,:) .* m(np,:) ./ m(1,:)));
hold off;
title(['Incandescence at ', num2str(np), ' ns v. fluence']);

figure(5);
cmap_sweep(nf, cm);
plot(t, X);
title('Annealed fraction');

figure(6);
plot(F0_vec, Ct(prompt, :), 'r');
hold on;
plot(F0_vec, Cinf(prompt, :), 'k');
hold off
ylim([0, 1.1]);
xlim([0, max(F0_vec)]);

