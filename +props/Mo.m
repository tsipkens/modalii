
% MO Define properties for molybdenum particles
% Author: Timothy Sipkens, 2019-11-03
%=========================================================================%

function prop = Mo(prop,opts)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = struct(); end

if ~exist('opts','var'); opts = struct(); end
if ~isfield(opts,'rho'); opts.rho = 'default'; end
if ~isfield(opts,'cp'); opts.cp = 'default'; end
if ~isfield(opts,'hv'); opts.hv = 'default'; end
if ~isfield(opts,'pv'); opts.pv = 'default'; end
if ~isfield(opts,'Em'); opts.Em = 'default'; end
%-------------------------------------------------------------------------%


%-- Sensible energy properties -------------------------------------------%
prop.phi = prop.h*prop.c/prop.kb;

prop.M = 0.09594;
prop.Tm = 2896; % (Ida and Guthrie, book)

switch opts.rho
    case {'default','Paradis'}
        prop.Arho = 1;
        prop.Brho = 1;
        prop.rho = @(T) 10789.^(4/3)./(props.iif(T>=prop.Tm,...
            (prop.Arho.*9100-prop.Brho.*0.6.*(T-prop.Tm)),...
            (prop.Arho.*9490-prop.Brho.*0.5.*(T-prop.Tm)))).^(1/3); % (Paradis,2002)
    case 'constant'
        prop.rho = @(T) 11500;
end

switch opts.cp
    case {'default','Paradis'}
        prop.Ccp = 1;
        prop.Dcp = 1;
        prop.Ecp = 1;
        prop.cp = @(T) props.iif(T>=prop.Tm,...
            (prop.Ccp.*34.2+prop.Dcp.*0.00113.*(T-prop.Tm))/prop.M,...
            (prop.Ccp.*44.0162+prop.Dcp.*0.0166.*(T-prop.Tm)+prop.Ecp.*5.5878e-07.*(T-prop.Tm).^2)./prop.M)...
            ; % (Paradis,2002)
        % prop.Ccp.*(151.78+0.00565.*(T-prop.Tm)).*(0.29+0.0000989.*(T-prop.Tm))./prop.M)
    case 'constant'
        prop.cp = @(T) 360;
end


%-- Conduction properties ------------------------------------------------%
prop.alpha = 0.15;
prop.ct = @(prop) sqrt(8 * prop.kb * prop.Tg / (pi * prop.mg));


%-- Evaporation properties -----------------------------------------------%
prop.mv = prop.M.*1.660538782e-24;							
prop.Rs = prop.R./prop.M;

prop.Tb = 4912; % (Ida and Guthrie, book)
prop.hvb = 590e3/1e6/prop.M; % (Ida and Guthrie, book)
switch opts.hv
    case {'default','Watson'}
        prop.n = 0.38; % Watson
        prop.Tcr = 14588; % (Young and Alder)
        prop.hvA = @()(prop.hvb*1e6)./...
            ((1-prop.Tb/prop.Tcr).^0.38); % Watson eqn. constant
        prop.hv = @(T) props.eq_watson(prop, T);
    case {'constant'}
        prop.hv = @(T) prop.hvb*1e6;
end

prop.Pref = 101325; % atmospheric boiling point used
prop.C = log(prop.Pref) + (prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
switch opts.pv
    case {'default','Kelvin-CC'}
        prop.gamma = @(dp,T) 2.11; % (currently unknown)
        prop.pv = @(T,dp,hv) props.eq_kelvin(prop, T,dp,hv);
    case {'Tolman-CC'}
        ... % enter additional parameters
        prop.gamma = @(dp,T) props.eq_tolman(prop, dp,T);
        prop.pv = @(T,dp,hv) props.eq_kelvin(prop, T,dp,hv);
    case {'CC'}
        prop.pv = @(T,dp,hv) props.eq_claus_clap(prop, T,dp,hv);
    case {'Antoine'}
        % will enter additional parameters
end


%-- Optical properties ---------------------------------------------------%
prop.CEmr = 1;
switch opts.Em
    case {'default','Barnes'}
        prop.Em_data = getfield(load('+props/Em_Mo_Barnes.mat'),'Em_data'); % (Barnes) *check
        prop.Em_gi = griddedInterpolant(prop.Em_data(:,1),prop.Em_data(:,2),'pchip');
        prop.Em = @(l,dp) prop.Em_gi(l);
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser);
    case {'Mie'}
        prop.Em_gi = PlasmaMie.Mo_Mie_inter;
        prop.Em = @(l,dp) prop.Em_gi(l,dp);
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1,dp)./prop.Em(l2,dp);
        prop.Eml = @(dp) prop.Em_gi(prop.l_laser,dp);
end


prop.opts = opts;

end

