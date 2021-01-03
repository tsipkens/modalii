
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

switch htmodel.opts.cond
    case 'free-molecular'
        q = q_fm(prop, T, dp, prop.Tg);  % compute free-molecular conduction
        
    case 'continuum'
        q = q_cont(prop, T, dp, prop.Tg);  % compute continuum regime conduction
        
    case {'transition', 'fuchs'}
        res = @(Tdelta) q_fm(prop, T, dp, Tdelta) - ...
            qcont(prop, Tdelta, dp + 2*get_mfp(prop), prop.Tg);
        
        Tdelta = fzero(@(T) res(T), [prop.Tg, T]);
        q = q_fm(prop, T, dp, Tdelta);

end

end



%== Q_FM =================================================================%
%   Free molecular conduction, using the thermal accommodation coefficient
%   and gamma2 (related to degrees of freedom) as a function of temperature.
function [q] = q_fm(prop, T, dp, Tg)

prop.alpha = min(max(prop.alpha, 0), 1); % added to force constraints
q = ((prop.alpha * prop.Pg * prop.ct(prop) * pi .* ...
    (dp .^ 2) ./ (8 * Tg)) .* ...
    prop.gamma2(T) .* (T - Tg));

end



%== Q_CONT ===============================================================%
%   Continuum regime condutions, using Eq. (6) in Daun and Hubermann. 
%   REQUIRES prop.k (conductivity) to be specified for the gas.
function [q] = q_cont(prop, T, dp, Tg)

q = 2 .* pi .* dp .* integral(@(T) prop.k(T), Tg, T);

end



%== GET_MFP ==============================================================%
%   Returns the Maxwell mean free path of the gas in nm based 
%   on viscosity (Eq. 16 in %Huberman and Daun) Pg is pressure in Pa, 
%   Tg is temperature in K, mg is molecular mass in kg.
%   REQUIRES prop.mu (dynamic viscosity) to be specified for the gas.
function lambda = get_mfp(prop)

rho = prop.mg * prop.Pg / (prop.kB * prop.Tg);

lambda = prop.mu / rho / sqrt(2*prop.kB .* prop.Tg / (pi * prop.mg));

end



