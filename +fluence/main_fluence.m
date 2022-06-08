
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

% prop = props.x_apb17_Fe;
prop = props.x_ldf;
prop = props.Ar(prop);
prop.Ti = 298;

% prop = props.Fe(prop, opts);
% prop = props.Ge(prop, opts); prop.Ti = 1675;
prop = props.C(prop, opts); prop.Ti = 1675;

prop.F0 = 0.15; % in [J/cm2]
prop.Tg = prop.Ti;
prop.sigma = 0;  % 0.1

% Define models and their parameterizations.
opts.abs = 'include';  % set opts to include absorption
x_fields = {'dp0', 'F0', 'CEmr'};  % set models to take only diameter as inputs
htmodel = HTModel(prop, x_fields, t, opts);
smodel = SModel(prop, x_fields, t, l);
smodel.htmodel = htmodel;
toc;
disp('Completed setup.');
disp(' ');


nf = 60;
F0_vec = linspace(0.0005, 2, nf);
dp = 30;
prop.dp0 = dp;

T = [];  J1 = [];  J2 = [];
disp('Computing temperature decays:');
tools.textbar([0, nf]);
for ii=1:length(F0_vec)
    [T(:,ii), m(:,ii)] = htmodel.evaluate([dp, F0_vec(ii), 1]);

    Jt = smodel.evaluateF([dp, F0_vec(ii), 1]);
    J1(:,ii) = Jt(:,1,1);  % choose first wavlength
    J2(:,ii) = Jt(:,1,1);  % choose second wavlength

    tools.textbar([ii, nf]);
end


%%
%{
figure(1);
cmap_sweep(nf, flipud(internet));
plot(t, T);
%}

propt = prop;
propt.hvb = prop.hvb / prop.M;  % convert from molar
Tlow1 = 10000 * 6 * pi * prop.Eml(dp) / ...
    (prop.l_laser * 1e-9 * prop.rho(prop.Tg) * prop.cp(prop.Tg)) ...
    .* F0_vec + prop.Tg;
[Tref, Fref] = fluence.calcTransition(propt);  % requires lambertw function
[Tfun, Tlow, Thigh] = fluence.getPeak(propt, [], [], Tref, Fref);

figure(2);
plot(F0_vec, max(T));
hold on;
plot(F0_vec, Tlow1, 'r');
plot(F0_vec, Tlow(F0_vec));
plot(F0_vec, Thigh(F0_vec));
plot(F0_vec, Tfun(F0_vec), 'k', 'LineWidth', 1.2); % plot overall fluence curve
hold off;
ylim([prop.Tg, 5500]);

%{
figure(3);
plot(F0_vec, max(J1));
hold on;
plot(F0_vec, max(J1 .* m ./ m(1,:)));
hold off;

np = 400;
figure(4);
plot(F0_vec, J1(np,:));
hold on;
plot(F0_vec, (J1(np,:) .* m(np,:) ./ m(1,:)));
hold off;
%}







