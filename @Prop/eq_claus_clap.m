
function [pv] = eq_claus_clap(prop,T,dp,hv)

pv = exp(prop.C-prop.hvb*1e6./prop.Rs./T); % Evaluate the Clausius-Clapeyron Eqn.

end