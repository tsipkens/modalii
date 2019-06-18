function [dXdt] = dXdt_fun(htmodel,q_ann,T,dp,X)
%DXDT Gets second output from Qann function.

[~, dXdt] = q_ann(T,dp,X);

end

