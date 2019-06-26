
% Q_COND Rate of free-molecular conduction energy loss from the nanoparticle.
% Author: Timothy Sipkens, 2018-12-17
%=========================================================================%

function [q] = q_cond(htmodel,T,dp)
%-------------------------------------------------------------------------%
% Inputs:
%   T       Vector of nanoparticle temperature, [K]
%   dp      Nanoparticle diameter, [nm]
%
% Outputs:
%   q       Rate of conductive losses, [W]
%-------------------------------------------------------------------------%

dp = dp.*1e-9; % convert to meters so everything is in SI units
prop = htmodel.prop;
prop.alpha = min(max(prop.alpha,0),1); % added to force constraints
q = ((prop.alpha*prop.Pg*prop.ct()*pi.*(dp.^2)./(8*prop.Tg)).*...
    prop.gamma2(T).*(T-prop.Tg));

end

