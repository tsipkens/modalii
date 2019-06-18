function [] = C(prop)

% Sensible energy properties **********************************************
prop.phi = prop.h*prop.c/prop.kb;
switch prop.opts.M % Molar mass in kg/mol
    case {'default','Kock','C3'}
        prop.M = 3*0.01201;
    case {'Liu','C2'}
        prop.M = 2*0.01201;
end

switch prop.opts.rho % Density in kg/m^3
    case {'default','Michelsen'}
        prop.Arho = 1;
        prop.Brho = 1;
        prop.rho = @(T) (prop.Arho.*2.303-prop.Brho.*7.3106e-5.*T).*1000; % Michelsen
    case 'Liu'
        prop.Arho = 1;
        prop.rho = @(T) prop.Arho.*1.9*1000.*ones(size(T));
    case {'Kock','constant'}
        prop.Arho = 1;
        prop.rho = @(T) prop.Arho.*1.86*1000;
end

switch prop.opts.cp % Specific heat in J/(kg K)
    case {'default','Michelsen'}
        prop.Ccp = 1;
        prop.cp = @(T) prop.Ccp.*(prop.R/0.01201).*(1.115.*(597./T).^2.*exp(597./T).*(exp(597./T)-1).^-2+...
            1.789.*(1739./T).^2.*exp(1739./T).*(exp(1739./T)-1).^-2+T./8260)...
            ; % Michelsen
    case 'Liu'
        Tt = 1200; % Transition temperature, Liu
        prop.Ccp = 1;
        prop.cp = @(T) prop.Ccp.*1000.*prop.iif(T<Tt,...
            -9.7768e-4+(2.7943e-4*T)+(1.4554e-5*T^2)-(3.4432e-8*T^3)+...
                (3.6700e-11*T^4)-(1.9485e-14*T^5)+(4.1802e-18*T^6),...
            2.9497e-1+(2.9614e-3*T)-(2.1232e-6*T^2)+(8.1901e-10*T^3)-...
                (1.7516e-13*T^4)+(1.9628e-17*T^5)-(8.9817e-22*T^6));
    case 'constant'
        prop.Ccp = 1;
        prop.cp = @(T) prop.Ccp.*1900;
end

% Conduction properties ***************************************************
prop.alpha = 0.37; % Liu
prop.ct = @()sqrt(8*prop.kb*prop.Tg/(pi*prop.mg));

% Evaporation properties **************************************************
switch prop.opts.mv % Molar mass in kg/mol
    case {'Liu'}
        prop.mv = @(T) (17.179+6.8654e-4.*T +2.9962e-6.*T.^2-...
            8.5954e-10.*T.^3+1.0486e-13.*T.^4)./1000.*1.660538782e-24;
    case {'piecewise-constant'}
        trans_fun = @(T,a,b) 1./(1+exp(-1./b.*(T-a)));
            % a sigmoid function is used to transition
        prop.mv = @(T) (12+trans_fun(T,1975,250).*24+...
            trans_fun(T,4200,250).*12)./1000.*1.660538782e-24;
    otherwise % 'C3', 'default', 'Kock'
        prop.mv = 3*0.01201.*1.660538782e-24;
end
						
prop.Rs = prop.R./prop.M;

prop.Tb = 4136.78; % Michelsen, C3 only
switch prop.opts.hv
    case {'default','Michelsen'}
        prop.hv = @(T) (8.443e5-26.921.*T)./prop.mv; % Michelsen, C3 only
    case 'Liu'
        prop.hv = @(T) (2.05398e5+(7.3660e2.*T)-(0.40713.*T.^2)+...
            (1.1992e-4.*T.^3)-(1.7946e-8.*T.^4)+...
            (1.0717e-12.*T.^5))./(prop.mv(T)./1.660538782e-24);
    case {'constant','Kock'}
        prop.hv = @(T) 7.9078e5.*ones(size(T))./prop.M; % Kock
    case {'Watson'}
        hfus = 29.3*3.4838e+05; % heat of fusion, Leider, 1973
        Tref = 4000;%5500; % reference point, Leider, 1973
        Nv = @(T) prop.mv(T)./prop.mv(Tref); % correct for J.kg conversion @ 5500K
        hvref = 11.5475;%7.8569; % reference point at T = 5500 K, Leider, 1973, J/kg @ 5500 K
        prop.Tcr = 6810; % Leider, 1973
        prop.n = 0.38;
        prop.hv = @(T) (((hvref*1e6).*...
            (((1-T./prop.Tcr)./(1-Tref/prop.Tcr)).^prop.n)))./Nv(T)+...
            hfus.*((T)<=4765); % 4765 K: triple point estimate, Leider, 1973
    case {'Roman'} % only valid for vapor
        hfus = 29.3*3.4838e+05; % heat of fusion, Leider, 1973, 30 kcal/gatom
        Tref = 4000;%5500; % reference point, Leider, 1973
        Nv = @(T) prop.mv(T)./prop.mv(Tref); % correct for J.kg conversion @ 4000K
        hvref = 11.5475;%7.8569; % reference point at T = 5500 K, Leider, 1973
        prop.Tcr = 6810; % Leider, 1973
        prop.n = 0.38; % Watson/Roman
        prop.beta = 0.371;
        prop.hv = @(T) (((hvref*1e6).*exp((prop.n-prop.beta).*...
            ((T-Tref)./(prop.Tcr-Tref))).*...
            (((1-T./prop.Tcr)./(1-Tref/prop.Tcr)).^prop.n)))./Nv(T)+...
            hfus.*((T)<=4765); % 4765 K: triple point estimate, Leider, 1973
end

switch prop.opts.pv % Vapor pressure in Pa, hvb in MJ to aid in stability
    case {'default','Michelsen-CC'}
        prop.hvb = prop.hv(prop.Tb)/1e6;
        prop.Pref = 101325; % atmospheric boiling point used
        prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
        prop.gamma = @(dp,T) 0.18; % (Shih, 2013)
        prop.pv = @prop.kelvinEqn;
    case 'Sipkens'
        prop.Tb = 4000;
        factor = ((((hvref*1e6).*exp((prop.n-prop.beta).*...
            ((prop.Tb-Tref)./(prop.Tcr-Tref))).*...
            (((1-prop.Tb./prop.Tcr)./(1-Tref/prop.Tcr)).^prop.n)))*...
            prop.mv(prop.Tb)/1.660538782e-24/(3*0.01201)+...
            hfus)/prop.hv(prop.Tb); % extra factor would correct for ~(Nv=3) at this point
        prop.hvb = prop.hv(prop.Tb)/1e6*1.03;%factor;
        prop.Pref = 101325*1.68;%103.4; % Leider et al., Table 3
        prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn.
        Tref2 = 4765;
        Pref2 = exp(prop.C-prop.hvb*1e6./prop.Rs./4765);
        hvb2 = prop.hvb-(29.3*3.4838e+5/1e6);
        C2 = log(Pref2)+(hvb2*1e6)./prop.Rs./Tref2; % Constant for C-C Eqn. 
        prop.gamma = @(dp,T) 0.18; % (Shih, 2013)
        pv0 = @(T) exp(prop.C-prop.hvb*1e6./prop.Rs./T).*(T<=4765)+... % 4765 K: triple point estimate, Leider, 1973
            exp(C2-hvb2*1e6./prop.Rs./T).*(T>4765);%exp(prop.C-prop.hvb*1e6./prop.Rs./T)*(T>4765); % Clausius-Clapeyron equation
        prop.pv = @(T,dp,hv) pv0(T).*...
            exp((4*prop.gamma(dp,T))./((dp).*prop.rho(T).*prop.Rs.*T)); % Kelvin equation
    case 'Liu'
        prop.hvb = prop.hv(prop.Tb)/1e6;
        prop.pv = @(T,dp,hv) 101325.*exp(-122.96+(9.0558e-2.*T)-(2.7637e-5.*T.^2)+...
            (4.1754e-9.*T.^3)-(2.4875e-13.*T.^4)); % Liu
    case {'Kock','constant'}
        prop.hvb = prop.hv(prop.Tb)/1e6;
        prop.Pref = 61.5;
        prop.Tb = 3000; % artificial, reference temperature
        prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
        prop.pv = @prop.clausClap;
    case 'Tolman-CC'
        ... % enter additional parameters
        prop.gamma = @prop.tolmanEqn;
        prop.pv = @prop.kelvinEqn;
    case 'Antoine'
        ... % enter additional parameters
end

% Optical properties ******************************************************
switch prop.opts.Em
    case {'default','constant'}
        prop.Em = @(l,dp) 0.34.*ones(1,length(l)); % Michelsen
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr;
        prop.Eml = @(dp) 0.34;
    case 'Liu'
        prop.Em = @(l,dp) 0.4.*ones(1,length(l)); % Liu, constant
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr;
        prop.Eml = @(dp) 0.4;
end

end

