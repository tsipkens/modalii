
% Q_RAD  Rate of radiative energy loss from the nanoparticle.
%        Includes a sub-function that evaluates Planck's law.
% AUTHOR: Timothy Sipkens, 2018-12-20
% NOTE: May include extrapolation for E(m) with restriciton that E(m) >= 0.
% 
% INPUTS:
%   T       Vector of nanoparticle temperature, [K]
%   dp      Nanoparticle diameter, scalar value [nm]
%
% OUTPUTS:
%   q       Rate of conductive losses, [W]
%   rad     Wavelength dependent spectrums, Cabs*Ibl
%   l       Wavelngths corresponding to rad, used in integration, [nm]
%=========================================================================%

function [q,rad,l] = q_rad(htmodel,T,dp)

dp = dp * 1e-9;
prop = htmodel.prop;

nl = 1000;
l_peak = 2.9e6 ./ T; % output in nm
l_min = 1 / 5 * min(l_peak);
l_max = min(2.5 .* max(l_peak), 5000);
    % caps maximum wavelength to 5000 nm preventing dl from being too large
dl = (l_max-l_min) / (nl-1);
l = (l_min:dl:l_max)'; % wavelengths to integrate over

Cabs = pi ^ 2 * (dp) ^ 3 ./ (l .* 1e-9) .* max(prop.Em(l, dp), 0); % absorption cross section
rad = bsxfun(@times, Cabs, blackbody(T, l)); % Cabs*Ibl

q = nansum(bsxfun(@times, rad, dl .* 1e-9)); % integrate over wavelength

end


%== BLACKBODY ============================================================%
% Sub-function that evaluates Planck's law for a given 
% temperature and wavelength.
% 
% INPUTS:
%   l     Wavelength, [nm]
%   T     Temperature, [K]
%-------------------------------------------------------------------------%
function out = blackbody(T,l)

p0 = (l .* 1e-9) * T;
p1 = 1 ./ (exp(0.0143877696 ./ (p0))); % Number is phi=h*c/kB
    % Wein's approximation used in forward model for scaling factor paper
p2 = 1.191042868e-16 ./ (l .* 1e-9) .^ 5; % Number is 2*h*c^2
out = bsxfun(@times, p1, p2);

end




