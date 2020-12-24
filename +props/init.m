
% INIT  Initialize a property structure.
% AUTHOR: Timothy Sipkens, 2020-12-24
%=========================================================================%

function [prop] = init()

prop.opts = struct();     % intialize options associated with properties

prop.h = 6.62606957e-34;  % Planck's constant [m^2.kg/s]
prop.c = 2.99792458e8;    % Speed of light in a vacuum [m/s]
prop.kb = 1.3806488e-23;  % Boltzmann constant [m^2.kg/s^2/K]
prop.R = 8.3144621;       % Universal gas constant [J/mol/K]

prop.phi = 0.0143877696;  % Constant for Planck's law, phi = h*c/kb

end

