
% J_EVAP    Simple bridging function to evaluate evaporation rate using q_evap function.
% Author:   Timothy Sipkens
%=========================================================================%

function [J] = J_evap(htmodel,T,dp)

[~,J] = q_evap(htmodel,T,dp);

end

