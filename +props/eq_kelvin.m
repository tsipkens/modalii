
function [pv] = eq_kelvin(prop, T, dp, hv)

pv0 = props.eq_claus_clap(prop, T, dp, hv); % Clausius-Clapeyron equation
pv = pv0 .* exp((4 * prop.gamma(dp, T, prop)) ./ ...
    ((dp) .* prop.rho(T,prop) .* prop.Rs .* T)); % Evaluate the Kelvin Eqn.

end

