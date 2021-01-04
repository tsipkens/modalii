
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
%   Kn      Knudsen number, []. (Optional.)
%=========================================================================%

function [q, Kn] = q_cond(htmodel, prop, T, dp, opts_cond)

% By default, use the opts from this instance of htmodel.
if ~exist('opts_cond', 'var'); opts_cond = []; end
if isempty(opts_cond); opts_cond = htmodel.opts.cond; end

% Convert to meters so everything is in SI units.
dp = dp .* 1e-9;

switch opts_cond
    case 'free-molecular'
        q = q_fm(prop, T, dp, prop.Tg);  % compute free-molecular conduction
        
    case 'continuum'
        q = q_cont(prop, T, dp, prop.Tg);  % compute continuum regime conduction
        
    case {'transition', 'fuchs'}
        % Define residual function, where
        % T_delta is an intermediate temperature at 
        % free molecular-continuum boundary (dp + 2*MFP).
        q = [];
        for ii=1:length(dp)
            res = @(T_delta) q_fm(prop, T, dp(ii), T_delta) - ...
                q_cont(prop, T_delta, dp(ii) + 2 * get_mfp(prop, T), prop.Tg);

            Tdelta = fzero(@(T) res(T), [prop.Tg, T]);
            q(ii) = q_fm(prop, T, dp(ii), Tdelta);
        end
        
end

% If necessary, also output Knudsen number.
% Requires that prop.mu (dynamic viscosity) to be specified for the gas.
if nargout>1
    Kn = get_mfp(prop, T) ./ (dp ./ 2);
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
%   REQUIRES: prop.k (conductivity) to be specified for the gas.
function [q] = q_cont(prop, T, dp, Tg)

q = 2 .* pi .* dp .* integral(@(T) prop.k(T), Tg, T);

end



%== GET_MFP ==============================================================%
%   Returns the Maxwell mean free path of the gas in nm based 
%   on viscosity (Eq. 16 in %Huberman and Daun) Pg is pressure in Pa, 
%   Tg is temperature in K, mg is molecular mass in kg.
%   REQUIRES: prop.mu (dynamic viscosity) to be specified for the gas.
function lambda = get_mfp(prop, T)

rho = prop.mg * prop.Pg ./ (prop.kb .* prop.Tg);
lambda = prop.mu(T) / rho / sqrt(2*prop.kb .* prop.Tg / (pi * prop.mg));

end



