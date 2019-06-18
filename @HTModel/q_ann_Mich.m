function [q,dXdt] = q_ann_Mich(htmodel,T,dp,X)
% Q_ANN_MICH Rate of annealing from Michelsen, 2003. 

dp = dp.*1e-9; % convert to meters so everything is in SI units
prop = htmodel.prop;
Na = 6.0221409e23; % Avagadro's number
Np = dp.^3.*pi.*prop.rho(T)./6./prop.M.*Na; % number of atoms in nanoparticle
Xd = 0.01; % initial defect density

A_dis = 1e18;
E_dis = 9.6e5;
k_dis = A_dis.*exp(-E_dis./prop.R./T); % disociation
A_int = 1e8;
E_int = 8.3e4;
k_int = A_int.*exp(-E_int./prop.R./T); % interstitial movement
A_vac = 1.5e17;
E_vac = 6.7e5;
k_vac = A_vac.*exp(-E_vac./prop.R./T); % vacancy movement

DH_int = -1.9e4; % interstitial movement
DH_vac = -1.4e5; % vacancy movement
Nd = (1-X).*(Xd.*Np); % number of defects

q = -(DH_int.*k_int+DH_vac.*k_vac).*Nd./Na;
dXdt = -1./(Xd.*Np).*(X.*Np./2.*k_dis-(k_int+k_vac).*Nd);

end



