
clear;
clc;
close all;

addpath cmap;
tic;

t = -0:0.1:100; % time, laser pulse centered at t = 0
l = [442,716]; % measurement wavelengths

opts = [];
% opts.hv = 'constant';
opts.Em = 'default'; %'Krishnan'; %'Mie-Krishnan';

prop = props.x_apb17_Fe;
prop = props.Ar(prop);

prop = props.Fe(prop, opts);
% prop = props.Ag(prop, opts);
% prop = props.C(prop, opts);

prop.F0 = 0.15; % in [J/cm2]
prop.Ti = 298;
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
F0_vec = linspace(0.0005, 0.35, nf);

T = [];  J1 = [];  J2 = [];
disp('Computing temperature decays:');
tools.textbar([0, nf]);
for ii=1:length(F0_vec)
    [T(:,ii), m(:,ii)] = htmodel.evaluate([30, F0_vec(ii), 1]);

    Jt = smodel.evaluateF([30, F0_vec(ii), 1]);
    J1(:,ii) = Jt(:,1,1);  % choose first wavlength
    J2(:,ii) = Jt(:,1,1);  % choose second wavlength

    tools.textbar([ii, nf]);
end


%%
figure(1);
cmap_sweep(nf, flipud(internet));
plot(t, T);

figure(2);
plot(F0_vec, max(T));

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




% [Fref, Tref] = fluence.calc_transition(prop);  % requires lambertw function
