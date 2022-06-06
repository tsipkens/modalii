
% DE_SOLVE  Solves ODEs for temperature, mass, etc. over time. 
% Author:   Timothy Sipkens, 2018-11-28
% 
% INPUTS:
%   prop    Heat transfer model property structure (used instead of htmodel.prop)
%   dp0     Vector of nanoparticle diameters for which ODEs are to be solved, [nm]
%
% OUTPUTS:
%   Tout    Time-resolved temperature with a single column per input diameter, [K]
%   dpo     Time-resolved nanoparticle diameter, same format as above, [nm]
%   mpo     Time-resolved nanoparticle mass, same format as above, [fraction]
%   Xo      Time-resolved anneealed fraction, same format as above, [fraction]

function [Tout, dpo, mpo, Xo] = de_solve(htmodel, prop, dp0)

Nd = length(dp0); % number of size classes to consider in solver

htmodel = htmodel.de_build;


%-- Initial conditions ---------------------------------------------------%
Ti = prop.Ti .* ones(Nd,1); % initial temperature, [K]

%-- Initial mass ----%
if isempty(prop.rho0)
	prop.rho0 = prop.rho(prop.Tg);
end
mass_conv = 1e21; % added for stability, converts mass to attogram (ag)
mpi = (prop.rho0 .* (dp0.*1e-9).^3 .* (pi/6)) .* mass_conv; % inital mass, [ag]

Xi = 0.*ones(Nd,1); % initial annealed fraction, [fraction]
%-------------------------------------------------------------------------%


%-- Starting point exception ---------------------------------------------%
%   Note: This is included to allow fitting to only part of a signal, 
%   while still allowing for polydispersity effects (see Mo in Sipkens et al.
%   2017).
if htmodel.t(1)>0.1 % gives initial condition at t=0 rather than t=t(1)
    t = [0;htmodel.t];
    opts_tadd = 1;
else
    t = htmodel.t;
    opts_tadd = 0;
end
%-------------------------------------------------------------------------%


%-- Setup ODE ------------------------------------------------------------%
%   Note: Switch chooses whether to couple three ODEs to also solve to the
%   annealed fraction or to only conisder two couples ODEs.
switch htmodel.opts.ann % whether to solve for annealed fraction
    case {'include', 'Michelsen', 'Sipkens'}
        dydt = @(t, y) real( ...
            [htmodel.dTdt(t, y(1:Nd), abs(y(Nd+1:2*Nd))./mass_conv, y(2*Nd+1:3*Nd)) .* 1e-9;...
             htmodel.dmdt(t, y(1:Nd), abs(y(Nd+1:2*Nd))./mass_conv, y(2*Nd+1:3*Nd)) .* mass_conv .* 1e-9;...
             htmodel.dXdt(t, y(1:Nd), abs(y(Nd+1:2*Nd))./mass_conv, y(2*Nd+1:3*Nd)) .* 1e-9]);
        yi = [Ti; mpi; Xi];
    otherwise
        dydt = @(t, y) real( ...
            [htmodel.dTdt(t, y(1:Nd), abs(y(Nd+1:2*Nd))./mass_conv) .* 1e-9; ...
             htmodel.dmdt(t, y(1:Nd), abs(y(Nd+1:2*Nd))./mass_conv) .* mass_conv .* 1e-9]);
                % Reframe dydt for ode solver, added real function
                %   added abs(y(2)) to force positive mass
                %   .*1e-9 is to convert denominator of dydt to ns for solver 
        yi = [Ti; mpi];
end
%-------------------------------------------------------------------------%


%-- Solve ODE ------------------------------------------------------------%
%   Note: Two solver methods are available: the built-in Runge-Kutta solver
%   from MATLAB or a simple implementation of Euler's method. 
switch htmodel.opts.deMethod
    case {'default', 'ode23s'} % use MATLAB Runge-Kutta solver
        
        if strcmp(htmodel.opts.abs,'include')
                % limit step size to ensure solver sees absorption
            opts.MaxStep = (prop.tlp)/2;
        else
            opts = [];
        end 
        
        [~,yo] = ode23s(dydt,t,yi,opts); % primary solver
        
        if ~(length(yo(:,1))==length(t))
            yo = [yo;...
                bsxfun(@times,[prop.Tg,yo(end,2)],...
                ones(length(t)-length(yo(:,1)),2))];
            disp('WARNING: Error in length of vector output by ode solver.');
        end
        Tout = max(yo(:,1:Nd),prop.Tg); % removes some of the ringing for fast decays
        mpo = yo(:,Nd+1:2*Nd)./mass_conv;
        
        switch htmodel.opts.ann % consider percentage annealed variable
            case {'include', 'Michelsen', 'Sipkens'}
                Xo = yo(:, 2*Nd + 1:3*Nd);
            otherwise
                Xo = zeros(size(Tout));
        end
        
    case 'Euler' % implementation of simple Euler integration
        dt = 0.2;
        t_eval = t(1):dt:t(end);
        T_eval = zeros(length(t_eval),1);
        m_eval = zeros(length(t_eval),1);
        X_eval = zeros(length(t_eval),1);
        T_eval(1) = yi(1);
        m_eval(1) = yi(2);
        X_eval(1) = yi(3);
        for ii=2:length(t_eval)
            dydt_ii = dydt(t_eval(ii),[T_eval(ii-1),m_eval(ii-1),X_eval(ii-1)]);
            T_eval(ii) = T_eval(ii-1)+dydt_ii(1).*dt;
            m_eval(ii) = m_eval(ii-1)+dydt_ii(2).*dt;
            X_eval(ii) = X_eval(ii-1)+dydt_ii(3).*dt;
        end
        
        Tout = interp1(t_eval,T_eval,t);
        mpo = interp1(t_eval,m_eval,t) ./ mass_conv;
        Xo = interp1(t_eval,X_eval,t);
        
    otherwise
        disp('deMethod not available.');
        
end
%-------------------------------------------------------------------------%


%-- Post-process results -------------------------------------------------%
if opts_tadd==1 % remove added initial point that was added if t(1)>0
    Tout = Tout(2:end,:);
    mpo = mpo(2:end,:);
    Xo = Xo(2:end,:);
end

dpo = ((6.*mpo) ./ (prop.rho(Tout).*pi)) .^ (1/3) .* 1e9;
    % calculate diameter over time
mpo = mpo ./ mpo(1); % output relative change in particle mass over time
%-------------------------------------------------------------------------%

end

