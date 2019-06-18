function [J] = J_evap(htmodel,T,dp)
% J_EVAP Evaluates evaporation rate using q_evap function. 

[~,J] = q_evap(htmodel,T,dp);

end

