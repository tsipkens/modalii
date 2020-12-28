
% AG Define properties for silver particles
% Author: Timothy Sipkens, 2019-11-03
%=========================================================================%

function prop = Ag(prop,opts)

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

prop.M = 0.10787;							
prop.Tm = 1234.9; % (Ida and Guthrie, book)

switch opts.rho
    case {'default'}
        prop.Arho = 1;
        prop.Brho = 1;
        prop.rho = @(T) prop.Arho.*10465-prop.Brho.*0.9067.*T; % (Kirshenbaum, 1961), given in kg/(m3)
end

switch opts.cp
    case {'default','Paradis'}
        prop.Ccp = 1;
        prop.cp = @(T) prop.Ccp.*33.46./prop.M; % (Iada and Guthrie, book)
end


%-- Conduction properties ------------------------------------------------%
prop.alpha = 0.3;
prop.Tg = 298;
prop.Pg = 101325;
prop.ct = @(prop) sqrt(8 * prop.kb * prop.Tg / (pi * prop.mg));


%-- Evaporation properties -----------------------------------------------%
prop.mv = prop.M.*1.660538782e-24;							
prop.Rs = prop.R./prop.M;

prop.Tb = 2435; % (Ida and Guthrie, book)
prop.hvb = 253e3/prop.M/1e6; % (Ida and Guthrie, book)
switch opts.hv
    case {'default','Watson'}
        prop.Tcr = 6410; % (Young and Alder)
        prop.n = 0.38; % Watson equations
        prop.hvA = @()(prop.hvb*1e6)./...
            ((1-prop.Tb/prop.Tcr).^0.38); % Watson eqn. constant
        prop.hv = @(T) props.eq_watson(prop, T);
    
    case {'constant'}
        prop.hv = @(T) prop.hvb;
end

prop.Pref = 101325; % atmospheric boiling point used
prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
switch opts.pv
    case {'default','Kelvin-CC'}
        prop.gamma = @(dp,T) (1.0994-0.0002*T); % (currently unknown)
        prop.pv = @(T,dp,hv) props.eq_kelvin(prop, T, dp, hv);
    
    case {'Tolman-CC'}
        ... % enter additional parameters
        prop.gamma = @(dp,T) props.eq_tolman(prop, T, dp);
        prop.pv = @(T,dp,hv) props.eq_kelvin(prop, T, dp, hv);
    
    case {'CC'}
        prop.pv = @(T,dp,hv) props.eq_claus_clap(prop, T,dp,hv);
    
    case {'Antoine'}
        ... % enter additional parameters
end


%-- Optical properties ---------------------------------------------------%
switch opts.Em
    case {'default','Drude'}
        prop.Em_data = getfield(load('Em_Ag_Drude.mat'),'Em_data'); % (Drude)
        prop.Em_gi = griddedInterpolant(prop.Em_data(:,1),prop.Em_data(:,2),'pchip');
        prop.Em = @(l,dp) prop.Em_gi(l);
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
end


%-- Particle size and signal properties ----------------------------------%
prop.dp0 = 20;
prop.sigma = 0; % Default monodisperse


prop.opts = opts;

end

