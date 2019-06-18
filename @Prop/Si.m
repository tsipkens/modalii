function [] = Si(prop)

% Sensible energy properties **********************************************
prop.phi = prop.h*prop.c/prop.kb;

prop.M = 0.02808;							
prop.Tm = 1687; % (Desai, 1985)

switch prop.opts.rho
    case {'default','Rhim'}
        prop.Arho = 1;
        prop.Brho = 1;
        prop.Crho = 1;
        prop.rho = @(T)(prop.Arho.*2.58-prop.Brho.*1.59e-4.*(T-prop.Tm)-...
            prop.Crho.*1.15e-7.*(T-prop.Tm).^2).*1000;
        % prop.rho = @(T)prop.iif(T>=prop.Tm,...
        %     (prop.Arho.*2.58-prop.Brho.*1.59e-4.*(T-prop.Tm)-prop.Crho.*1.15e-7.*(T-prop.Tm).^2).*1000,...
        %     prop.Arho.*2329.*ones(length(T),1)); % (Rhim, 2000)
    case {'Assael'}
        prop.Arho = 1;
        prop.Brho = 1;
        prop.rho = @(T)prop.Arho.*2550-prop.Brho.*0.264.*(T-prop.Tm);
end

switch prop.opts.cp
    case {'default','Deasi'}
        prop.Ccp = 27.2/prop.M; % (Desai,1985), given in J/(kg K)
        prop.cp = @(T) prop.Ccp;
end

% Conduction properties ***************************************************
prop.alpha = 0.36;
prop.Tg = 298;
prop.Pg = 101325;
prop.ct = @()sqrt(8*prop.kb*prop.Tg/(pi*prop.mg));

% Evaporation properties **************************************************
prop.mv = prop.M.*1.660538782e-24;							
prop.Rs = prop.R./prop.M;

prop.Tb = 3538; % (CRC Handbook)
prop.hvb = 395e3/prop.M/1e6; % (Ida and Guthrie, book)
switch prop.opts.hv
    case {'default','Watson'}
        prop.Tcr = 5193; % (Ullmann's Encyclopedia of Industrial Chemistry)
        prop.n = 0.38; % Watson
        prop.hvA = @()(prop.hvb*1e6)./...
            ((1-prop.Tb/prop.Tcr).^0.38); % Watson eqn. constant
        prop.hv = @prop.watsonEqn;
    case {'constant'}
        prop.hv = @(T) prop.hvb;
end

prop.Pref = 101325; % atmospheric boiling point used
prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
switch prop.opts.pv
    case {'default','Kelvin-CC'}
        prop.gamma = @(dp,T) 0.765-0.016e-3.*(T-prop.Tm); % Rhim, 2000
        prop.pv = @prop.kelvinEqn;
    case {'Tolman-CC'}
        prop.delta = 0.22; % Toman length = atomic diameter
        prop.gamma0 = @(dp,T) 0.765-0.016e-3.*(T-prop.Tm); % Rhim, 2000
        prop.gamma = @prop.tolmanEqn;
        prop.pv = @prop.kelvinEqn;
    case {'CC'}
        prop.pv = @prop.clausClap;
    case {'Antoine'}
        ... % enter additional parameters
end

% Optical properties ******************************************************
switch prop.opts.Em
    case {'default','constant'}
        prop.CEmr = 1;
        prop.Em = @(l,dp) 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
end

% Particle size and signal properties *************************************
prop.C_dp = 100;
prop.dp0 = 20; % prop.alpha.*prop.C_dp
prop.sigma = 0; % Default monodisperse

end

