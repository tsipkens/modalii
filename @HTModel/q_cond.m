
% Q_COND Rate of free-molecular conduction energy loss from the nanoparticle.
% AUTHOR: Timothy Sipkens, 2018-12-17
% Contains code snippets from K. J. Daun.
% 
% INPUTS:
%   T       Vector of nanoparticle temperature, [K]
%   dp      Nanoparticle diameter, [nm]
%
% OUTPUTS:
%   q       Rate of conductive losses, [W]
%=========================================================================%

function [q] = q_cond(htmodel, prop, T, dp)

dp = dp .* 1e-9; % convert to meters so everything is in SI units
q = q_fm(prop, T, dp);  % compute free-molecular conduction

end



%== Q_FM =================================================================%
%   Free molecular conduction, using the thermal accommodation coefficient
%   and gamma2 (related to degrees of freedom) as a function of temperature.
function [q] = q_fm(prop, T, dp)

prop.alpha = min(max(prop.alpha, 0), 1); % added to force constraints
q = ((prop.alpha * prop.Pg * prop.ct(prop) * pi .* ...
    (dp .^ 2) ./ (8 * prop.Tg)) .* ...
    prop.gamma2(T) .* (T - prop.Tg));

end


%== Q_CONT ===============================================================%
%   Continuum regime condutions, using Eq. (6) in Daun and Hubermann. 
%   Requires prop.k to be specified for the given material.
%   NOTE: Currently not used.
function [q] = q_cont(prop, T, dp)

q = 2 .* pi .* dp .* integral(@(T) prop.k(T), prop.Tg, T);

end


