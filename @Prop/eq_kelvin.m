
function [pv] = eq_kelvin(prop, T, dp, hv)

pv0 = prop.eq_claus_clap(T,dp,hv); % Clausius-Clapeyron equation
pv = pv0.*exp((4*prop.gamma(dp,T))./((dp).*prop.rho(T).*prop.Rs.*T)); % Evaluate the Kelvin Eqn.

end

