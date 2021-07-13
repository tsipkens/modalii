
% SI  Define properties for silicon particles
%  
%  AUTHOR: Timothy Sipkens, 2019-11-03
%  
% _________________________________________________________________________

function prop = Si(prop, opts)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = props.init(); end

if ~exist('opts', 'var'); opts = struct(); end
if ~isfield(opts, 'rho'); opts.rho = 'default'; end
if ~isfield(opts, 'cp'); opts.cp = 'default'; end
if ~isfield(opts, 'hv'); opts.hv = 'default'; end
if ~isfield(opts, 'pv'); opts.pv = 'default'; end
if ~isfield(opts, 'Em'); opts.Em = 'default'; end
%-------------------------------------------------------------------------%


prop.matl = 'Si';


%-- Sensible energy properties -------------------------------------------%
prop.M = 0.02808;							
prop.Tm = 1687; % (Desai, 1985)

switch opts.rho
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

switch opts.cp
    case {'default','Deasi'}
        prop.Ccp = 27.2/prop.M; % (Desai,1985), given in J/(kg K)
        prop.cp = @(T) prop.Ccp;
end


%-- Conduction properties ------------------------------------------------%
prop.alpha = 0.36;
prop.Tg = 298;
prop.Pg = 101325;
prop.ct = @(prop) sqrt(8 * prop.kb * prop.Tg / (pi * prop.mg));


%-- Evaporation properties -----------------------------------------------%
prop.mv = prop.M.*1.660538782e-24;							
prop.Rs = prop.R./prop.M;

prop.Tb = 3538; % (CRC Handbook)
prop.hvb = 395e3/prop.M/1e6; % (Ida and Guthrie, book)
switch opts.hv
    case {'default','Watson'}
        prop.Tcr = 5193; % (Ullmann's Encyclopedia of Industrial Chemistry)
        prop.n = 0.38; % Watson
        prop.hvA = @()(prop.hvb*1e6)./...
            ((1-prop.Tb/prop.Tcr).^0.38); % Watson eqn. constant
        prop.hv = @(T) prop.eq_watson(T);
    case {'constant'}
        prop.hv = @(T) prop.hvb;
end

prop.Pref = 101325; % atmospheric boiling point used
prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
switch opts.pv
    case {'default','Kelvin-CC'}
        prop.gamma = @(dp,T) 0.765-0.016e-3.*(T-prop.Tm); % Rhim, 2000
        prop.pv = @(T,dp,hv) prop.eq_kelvin(T,dp,hv);
    case {'Tolman-CC'}
        prop.delta = 0.22; % Toman length = atomic diameter
        prop.gamma0 = @(dp,T) 0.765-0.016e-3.*(T-prop.Tm); % Rhim, 2000
        prop.gamma = @(dp,T) prop.eq_tolman(dp,T);
        prop.pv = @(T,dp,hv) prop.eq_kelvin(T,dp,hv);
    case {'CC'}
        prop.pv = @(T,dp,hv) prop.eq_claus_clap(T,dp,hv);
    case {'Antoine'}
        ... % enter additional parameters
end


%-- Optical properties ---------------------------------------------------%
switch opts.Em
    case {'default','constant'}
        prop.CEmr = 1;
        prop.Em = @(l,dp) 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
end


%-- Particle size and signal properties ----------------------------------%
prop.C_dp = 100;
prop.dp0 = 20; % prop.alpha.*prop.C_dp
prop.sigma = 0; % Default monodisperse


prop.opts = opts;

end

