function [q,h] = sig_q(htmodel,T,dp,opts_rad)
% SIG_Q Evaluates and plots the different heat transfer modes for the current model.
% Author: Timothy Sipkens, 2018-12-17
%
% Includes a sub-function that evaluates Planck's law.
% Note: May include extrapolation for E(m) with restriciton that E(m) >= 0.
%-------------------------------------------------------------------------%
% Inputs:
%   T           Vector of nanoparticle temperature, [K]
%   dp          Nanoparticle diameter, scalar value [nm]
%   opts_rad    Whether to plot the radiation mode, boolean     (Optional, default is true)
%
% Outputs:
%   q           Rate of losses, one column per mode [W]
%   h           Figure handle
%-------------------------------------------------------------------------%


if nargin<4
    opts_rad = true;
end

if opts_rad
    q = [q_cond(htmodel,T,dp)' q_evap(htmodel,T,dp)'-(J_evap(htmodel,T,dp).*T.*htmodel.prop.cp(T))' ...
        q_rad(htmodel,T,dp)'];
else
    q = [q_cond(htmodel,T,dp)' q_evap(htmodel,T,dp)'-(J_evap(htmodel,T,dp).*T.*htmodel.prop.cp(T))'];
end

%-- Plot modes -----------------------------------------------------------%
h = semilogy(T,q);

hold on; % plot vertical lines about transition points
fun1 = @(T)(q_cond(htmodel,T,dp)-q_evap(htmodel,T,dp)+(J_evap(htmodel,T,dp).*T.*htmodel.prop.cp(T)));
t1 = fzero(fun1,3000);
plot([t1,t1],ylim,'k--');

if opts_rad
    fun2 = @(T)(q_evap(htmodel,T,dp)-(J_evap(htmodel,T,dp).*T.*htmodel.prop.cp(T))-q_rad(htmodel,T,dp));
    t2 = fzero(fun2,3000);
    plot([t2,t2],ylim,'k--');
    legend('q_{cond}','q_{evap}','q_{rad}');
else
    legend('q_{cond}','q_{evap}');
end
hold off;

xlim([min(T),max(T)]); % so as to cut off vertical lines if out of scope of T



end

