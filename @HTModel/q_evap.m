
% Q_EVAP  Rate of evaporation or sublimation energy loss from the nanoparticle. 
% AUTHOR: Timothy Sipkens, 2018-12-17
% 
% INPUTS:
%   T       Vector of nanoparticle temperature, [K]
%   dp      Nanoparticle diameter, [nm]
%
% OUTPUTS:
%   q       Rate of evaporative or sublimative losses, [W]
%=========================================================================%

function [q,J,hv,pv] = q_evap(htmodel, prop, T, dp)

dp = dp .* 1e-9; % convert to meters so everything is in SI units
prop = htmodel.prop;

if isempty(prop.gamma)
	prop.gamma = @props.eq_tolman;
end

if isempty(prop.alpham)
	prop.alpham = @(T) 1;
end

hv = prop.hv(T);
pv = prop.pv(T,dp,hv);
if isa(prop.mv,'function_handle')
    mv = prop.mv(T);
else
    mv = prop.mv;
end

cv = sqrt(8 * prop.kb .* T ./ (pi * mv)); % Molecular speed of matl [m/s]
nv = prop.alpham(T) .* pv ./ (prop.kb .* T);  % Vapor number flux of matl

J = mv .* nv .* cv ./ 4 .* pi .* dp .^ 2;
q = hv .* J;

end

