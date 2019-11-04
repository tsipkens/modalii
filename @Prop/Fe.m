function prop = Fe(prop)

% Sensible energy properties **********************************************
prop.phi = prop.h*prop.c/prop.kb;

prop.M = 0.055847;
prop.Tm = 1811; % (Desai, 1986)
switch prop.opts.rho
    case {'default','Hixson'} % (Hixson, 1990), T > 2125 K
        prop.Arho = 8171;
        prop.Brho = -0.64985;
        prop.rho = @(T) (prop.Brho.*T+prop.Arho);
        
    case {'Basinski'} % (Basinksi,1955 )
        prop.rho_data = getfield(load('rho_Fe_Basinksi.mat'),'rho_data');
        prop.rho_gi = griddedInterpolant(prop.rho_data(:,1),prop.rho_data(:,2),'pchip');
        prop.rho = @(T) prop.rho_gi(T);
        
    case {'Mills'}
        prop.rho = @(T)...
            prop.iif(and(T>=0,T<293),7871.6.*(T./T),...
            prop.iif(and(T>=293,T<1184),7874./(1+3.*(14.5e-6).*(T-273-20)),...
            prop.iif(and(T>=1184,T<1667),7650-0.51.*(T-273-911),...
            prop.iif(and(T>=1667,T<1811),7355-0.42.*(T-273-1394),...
            7030-0.86.*(T-273-1538)))));
        
    case 'noSlope'
        prop.Arho = 6350; % 7874
        prop.Brho = 0;
        prop.rho = @(T) (prop.Brho.*T+prop.Arho);
        
    case 'constant'
        prop.Arho = 6350; % 7874
        prop.rho = @(T) prop.Arho;
        
end

switch prop.opts.cp
    case {'default','mixed'}
        prop.Ccp = 1;
        prop.Dcp = 1; % update in the future, originally zero for APB paper
        prop.cp = @(T) prop.Ccp.*(prop.iif(T>=prop.Tm,(46.632.*ones(1,length(T))),... % liquid iron, (Desai, 1986)
            prop.iif(T>=1667,(-12.38+3.161e-2.*T),... % delta iron, (Cezairliyan and McClure, 1974)
            (17.64+1.232e-2.*T)))./prop.M)... % gamma iron, (Cezairliyan and McClure, 1974)
            ; % given in J/(kg K)
        
    case 'Desai'
        prop.Ccp = 1;
        prop.Dcp = 1;
        prop.cp = @(T) (prop.iif(T>=prop.Tm,(prop.Ccp.*46.632.*ones(1,length(T))),... % liquid iron, (Desai, 1986)
            prop.iif(T>=1667,(prop.Ccp.*40.368+prop.Dcp.*3.2194e-2.*(T-1667)),... % delta iron, (Desai, 1986)
            (prop.Ccp.*33.803+prop.Dcp.*9.1605e-3.*(T-1181))))./prop.M)... % gamma iron, (Desai, 1986)
            ; % given in J/(kg K)
        
    case 'noSlope'
        prop.Ccp = 1; % 7874
        prop.Dcp = 0;
        prop.cp = @(T) prop.Ccp.*(prop.iif(T>=prop.Tm,(46.632.*ones(1,length(T))),... % liquid iron, (Desai, 1986)
            prop.iif(T>=1667,(-12.38+3.161e-2.*T),... % delta iron, (Cezairliyan and McClure, 1974)
            (17.64+1.232e-2.*T)))./prop.M)... % gamma iron, (Cezairliyan and McClure, 1974)
            ; % given in J/(kg K)
        
    case 'constant'
        prop.Ccp = 1;
        prop.cp = @(T) prop.Ccp.*46.632./prop.M; % given in J/(kg K)
end

% Conduction properties ***************************************************
prop.alpha = 0.23; % 0.23 - model selection						
prop.Tg = 298;					
prop.Pg = 101325;													
prop.ct = @()sqrt(8*prop.kb*prop.Tg/(pi*prop.mg));

% Evaporation properties **************************************************
prop.mv = prop.M.*1.660538782e-24;
prop.Rs = prop.R./prop.M;

prop.Tb = 3134; % (Dean, Lange's Handbook of Chemisty) *check
prop.hvb = 340e3/prop.M/1e6; % (Dean, Lange's Handbook of Chemisty) *check
switch prop.opts.hv
    case {'default','Watson'}
        prop.Tcr = 9340; % (Young and Alder) (ALT: Beutl et al.)
        prop.n = 0.38; % Watson
        prop.hvA = @()(prop.hvb*1e6)./...
            ((1-prop.Tb/prop.Tcr).^prop.n); % Watson eqn. constant
        prop.hv = @prop.watsonEqn;
        
    case {'Roman'}
        prop.Tcr = 9340; % (Young and Alder) (ALT: Beutl et al.)
        prop.n = 0.38; % Watson/Roman
        prop.beta = 0.371;
        prop.hv = @(T) (prop.hvb*1e6).*exp((prop.n-prop.beta).*...
            ((T-prop.Tb)./(prop.Tcr-prop.Tb))).*...
            (((prop.Tcr-T)./(prop.Tcr-prop.Tb)).^prop.n);
        
    case {'Meyra'}
        prop.Tcr = 9340; % (Young and Alder) (ALT: Beutl et al.)
        prop.n = 0.38; % Watson/Roman
        prop.hv = @(T) (prop.hvb*1e6).*...
            (((prop.Tcr-T)./(prop.Tcr-prop.Tb))...
            .^((prop.n^2).*((prop.Tcr-T)./(prop.Tcr-prop.Tb))+prop.n));
        
    case {'constant'}
        prop.hv = @(T) prop.hvb*1e6;
end

prop.gamma0 = 1.865; % (Keene et al, 1988)
prop.Pref = 101325; % atmospheric boiling point used
prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % constant for C-C Eqn. 
switch prop.opts.pv
    case {'default','Kelvin-CC'}
        % prop.gamma = @(dp,T)(prop.gamma0-0.35*1e-3*(T-1823)); % (Keene et al, 1988)
        prop.gamma = @(dp,T) prop.gamma0;
        prop.pv = @prop.kelvinEqn;
        
    case {'Tolman-CC'}
        ... % enter additional parameters
        % prop.gammaT = @(T)(prop.gamma0-0.35*1e-3*(T-1823)); % (Keene et al, 1988)
        prop.gammaT = @(T) prop.gamma0;
        prop.delta = 0.126; % atomic diameter
        prop.gamma = @prop.tolmanEqn;
        prop.pv = @prop.kelvinEqn;
        
    case {'CC'}
        prop.pv = @prop.clausClap;
        
    case {'CC-alt'}
        prop.C1 = prop.hvb*1e6./prop.Rs;
        pv0 = exp(prop.C-prop.C1./T);
        prop.pv = pv0.*exp((4*prop.gamma0)./...
            ((dp).*prop.rho(T).*prop.Rs.*T));
        
    case {'Antoine-alt'}
        prop.C1 = prop.hvb*1e6./prop.Rs;
        prop.C2 = 1;
        pv0 = exp(prop.C-prop.C1./(T+prop.C2));
        prop.pv = pv0.*exp((4*prop.gamma0)./...
            ((dp).*prop.rho(T).*prop.Rs.*T));
        
    case {'Rankine-Kirchoff-alt'}
        prop.C1 = prop.hvb*1e6./prop.Rs;
        prop.C2 = 1;
        pv0 = exp(prop.C-prop.C1./T+prop.C2.*log(T));
        prop.pv = pv0.*exp((4*prop.gamma0)./...
            ((dp).*prop.rho(T).*prop.Rs.*T));
        
    case {'Nernst-alt'}
        prop.C1 = prop.hvb*1e6./prop.Rs;
        prop.C2 = 1;
        prop.C3 = 1;
        pv0 = exp(prop.C-prop.C1./T+prop.C2.*log(T)+prop.C3.*T);
        prop.pv = pv0.*exp((4*prop.gamma0)./...
            ((dp).*prop.rho(T).*prop.Rs.*T));
        
    case {'CRC-alt'}
        prop.C1 = prop.hvb*1e6./prop.Rs;
        prop.C2 = 1;
        prop.C3 = 1;
        prop.C4 = 1;
        pv0 = exp(prop.C-prop.C1./T+prop.C2.*log(T)+prop.C3./(T.^3));
        prop.pv = pv0.*exp((4*prop.gamma0)./...
            ((dp).*prop.rho(T).*prop.Rs.*T));
        
    case {'Kelvin-Antoine'}
        prop.gamma = @(dp,T) prop.gamma0;
        pv0 = prop.Antoine(T,dp,hv); % Clausius-Clapeyron equation
        prop.pv = pv0.*exp((4*prop.gamma(dp,T))./((dp).*prop.rho(T).*prop.Rs.*T));
            % Evaluate the Kelvin Eqn.
end

% Optical  properties *****************************************************
switch prop.opts.Em
    case {'default','Emr1.1'}
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*1.1;
        prop.Em = @(l,dp) (l-716)/(442-716)*(prop.Emr(442,716)-1)+1;
        
    case {'Drude'}
        prop.omega_p = 6.78e17;
        prop.tau = 1.69e-19;
        % [Em_temp,n_temp,k_temp] = prop.drude(400:1100);
        prop.Em_data = getfield(load('Em_Fe_Drude.mat'),'Em_data');
        prop.Em_gi = griddedInterpolant(prop.Em_data(:,1),prop.Em_data(:,2),'pchip');
        prop.Em = @(l,dp) prop.Em_gi(l);
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1,dp)./prop.Em(l2,dp);
        
    case {'Krishnan'}
        prop.Em = @(l,dp) ones(size(dp))*...
            (polyval([3.9751e-13,-1.6904e-9,2.7217e-6,...
            -0.0020557,0.69807],l)); % quartic fit to Krishnan data
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1,dp)./prop.Em(l2,dp);
        
    case {'Shvarev'}
        prop.Em_data = getfield(load('Em_Fe_Shvarev.mat'),'Em_data');
        prop.Em_gi = griddedInterpolant(prop.Em_data(:,1),prop.Em_data(:,2),'pchip');
        prop.Em = @(l,dp) prop.Em_gi(l);
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1,dp)./prop.Em(l2,dp);
        
    case 'Mie-Krishnan'
        load('@Prop\Fe_opt_prop.mat');
        n_gi = griddedInterpolant(Krishnan.l,Krishnan.n,'pchip','linear');
        k_gi = griddedInterpolant(Krishnan.l,Krishnan.k,'pchip','linear');
        n = @(l) n_gi(l);
        k = @(l) k_gi(l);
        prop.Em = @(l,dp) prop.get_Mie_solution(n,k,l,dp);
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1,dp)./prop.Em(l2,dp);
        
end
% prop.Eml = @(dp) 0.07; % representative of (Sipkens et al., APB, 2017)
prop.Eml = @(dp) prop.Em(prop.l_laser,dp);

% Particle size and signal properties *************************************
prop.dp0 = 20;
prop.sigma = 0; % Default monodisperse

end

