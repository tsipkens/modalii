
% CONVERTPROP  Updates prop format for use with +fleunce functions.
%  
%  This function converts prop structures from those used in other
%  parts of the LII code developed by T. A. Sipkens.

function prop = convertProp(prop_old)

%-- Experimental parameters ----------------------------------------------%
prop.Tg = prop_old.Tg; % gas temperature, K
prop.dp = prop_old.dp0; % nanoparticle diameter, m
prop.tlp = prop_old.tlp; % FWHM of the laser pulse, s


%-- Spectroscopic properties ---------------------------------------------%
prop.Eml = prop_old.Eml; % absorption function @ laser wavelength
prop.l_laser = prop_old.l_laser; % laser wavelength, m


%-- Sensible energy properties -------------------------------------------%
prop.rho = prop_old.rho(prop.Tg); % density, kg/m^3
prop.cp = prop_old.cp(prop.Tg); % specific heat capacity, J/(kg*K)
% Note: these expressions are evaluated the gas temperature.


%-- Vaporization properties ----------------------------------------------%
prop.M = prop_old.M; % molar mass of the vapor, kg/mol
prop.mv = prop.M*1.660538782e-24; % mass of the vapor, kg
prop.kb = 1.38064852e-23; % Boltzmann constant, m^2*kg/(s^2*K)
prop.R = 8.3144598; % unviersal gas constant, kg*?m^2/(s^2*K*?mol)
prop.Rs = prop.R/prop.M; % specific gas constant, ?m^2/(s^2*K)
prop.beta = prop_old.alpham();
if isempty(prop.beta); prop.beta = 1; end  % if not available in prop_old
prop.hvb = prop_old.hvb;
    % latent heat of vaporization at boiling point, J/kg

% Clausius-Clapeyron (C-C) equation.
prop.Tb = prop_old.Tb; % reference temperature of C-C equation, K
prop.Pb = prop_old.Pref; % reference pressure for C-C equation, Pa
prop.A = exp(log(prop.Pb)+prop.hvb/prop.Rs/prop.Tb); 
    % constant for C-C equation, Pa
% prop.pv = @(T) prop.A.*exp(-prop.hvref/prop.Rs./T); % C-C equation, Pa

% Parameter C1, as defined in Sipkens and Daun (2017).
prop.C1 = prop.kb * pi / (9 * prop.Rs * prop.hvb * prop.mv) * ...
    (prop.rho * prop.cp / (prop.A * prop.beta)) ^ 2; 


end