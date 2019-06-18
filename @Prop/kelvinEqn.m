function [pv] = kelvinEqn(prop,T,dp,hv)

pv0 = prop.clausClap(T,dp,hv); % Clausius-Clapeyron equation
pv = pv0.*exp((4*prop.gamma(dp,T))./((dp).*prop.rho(T).*prop.Rs.*T)); % Evaluate the Kelvin Eqn.

end

