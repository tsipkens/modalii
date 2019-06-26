
% DXDT Bridging function to get second output from Qann function.
% Author: Timothy Sipkens
%=========================================================================%

function [dXdt] = dXdt_fun(htmodel,q_ann,T,dp,X)

[~, dXdt] = q_ann(T,dp,X);

end

