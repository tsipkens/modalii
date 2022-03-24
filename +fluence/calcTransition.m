
% CALC_TRANSITION  Calculates fluence regime transition point. 
%  
%  This function takes a property struct and outputs the reference or
%  transition fluence (Fref) and temperature (Tref). 
%  
%  See Sipkens and Daun (Opt. Express, 2017) for more information. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2017

function [Tref, Fref] = calcTransition(prop)

prop = fluence.convertProp(prop);

min_fun = @(Tref) norm(...
    2 * prop.hvb / prop.Rs+...
    Tref .* lambertw(-1, -prop.C1 .* (prop.dp .* (Tref - prop.Tg) ./ prop.tlp) .^2 )...
    ); % function to minimize to solve for Tref

T0 = 3000; % starting point for minimizer, K
Tref = fminsearch(min_fun,T0); % calculate reference temperature, K

Fref = prop.l_laser * prop.rho * prop.cp * (Tref - prop.Tg) / ...
    (6 * pi * prop.Eml) / 10000; % /10000 converts fluence to J/cm^2

end

