
% Ge Define properties for germanium particles
% Author: Timothy Sipkens, 2019-11-03
%=========================================================================%

function prop = Ge(prop,opts)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('prop','var'); prop = struct(); end

if ~exist('opts','var'); opts = struct(); end
if ~isfield(opts,'rho'); opts.rho = 'default'; end
if ~isfield(opts,'cp'); opts.cp = 'default'; end
if ~isfield(opts,'hv'); opts.hv = 'default'; end
if ~isfield(opts,'pv'); opts.pv = 'default'; end
if ~isfield(opts,'Em'); opts.Em = 'default'; end
%-------------------------------------------------------------------------%


prop.matl = 'Ge';


%-- Sensible energy properties ---------------------------------------------%
prop.phi = prop.h*prop.c/prop.kb;

prop.M = 0.07264;
prop.Tm = 1210; % (Ida and Guthrie, book)

switch opts.rho
    case {'default'}
        prop.Arho = 5490; % (Ida and Guthrie, book)
        prop.Brho = -0.49; % (Ida and Guthrie, book)
        prop.rho = @(T) prop.Brho.*(T-prop.Tm)+prop.Arho;
    case {'Jan','Menser'}
        prop.Arho = 6170;
        prop.Brho = -4.42e-4;
        prop.rho = @(T) (prop.Brho.*T+prop.Arho);
    case {'constant'}
        prop.Arho = 5490;
        prop.rho = @(T) prop.Arho;
end

switch opts.cp
    case {'default','mixed','constant'}
        prop.Ccp = 29.3./prop.M*1.2; % (Ida and Guthrie, book)
        prop.Dcp = 0;
        prop.cp = @(T) prop.Ccp+prop.Dcp.*T;
    case {'Jan'}
        prop.Ccp = 971.42; %971.42; 1.1657e3;
        prop.cp = @(T) prop.Ccp;
end


%-- Conduction properties ------------------------------------------------%
prop.alpha = 0.23;
prop.ct = @(prop) sqrt(8 * prop.kb * prop.Tg / (pi * prop.mg));


%-- Evaporation properties -----------------------------------------------%
prop.mv = prop.M.*1.660538782e-24;
prop.Rs = prop.R./prop.M;

prop.Tb = 3103; % (Ida and Guthrie, book)
prop.hvb = 333e3/prop.M/1e6; % (Ida and Guthrie, book)
switch opts.hv
    case {'default','Watson'}
        prop.Tcr = 9803; % (Young and Alder)
        prop.n = 0.38;
        prop.hvA = @()(prop.hvb*1e6)./...
            ((1-prop.Tb/prop.Tcr).^0.38); % Watson eqn. constant
        prop.hv = @prop.watson;
    case {'constant'}
        prop.hv = @(T) prop.hvb*1e6;
    case {'Jan'}
        prop.hvb = 5.4;%4.256177038777618;
        prop.hv = @(T) prop.hvb*1e6;
end

prop.Pref = 101325; % atmospheric boiling point used
prop.C = log(prop.Pref)+(prop.hvb*1e6)./prop.Rs./prop.Tb; % Constant for C-C Eqn. 
switch opts.pv
    case {'default','Kelvin-CC'}
        prop.gamma = @(dp,T) 607e-3+(T-prop.Tm).*0.14e-3; % (Ida and Guthrie, book)
        prop.pv = @(T,dp,hv) prop.eq_kelvin(T,dp,hv);
    case {'Tolman-CC'}
        ... % enter additional parameters
        prop.gamma = @(dp,T) prop.eq_tolman(dp,T);
        prop.pv = @(T,dp,hv) prop.eq_kelvin(T,dp,hv);
    case {'CC'}
        prop.pv = @(T,dp,hv) prop.eq_claus_clap(T,dp,hv);
    case {'Antoine','Jan'}
        prop.C1 = 16149;
        prop.C2 = 11.7288;
        prop.pv = @(T,dp,hv) prop.eq_antoine(T,dp,hv);
end


%-- Optical properties ---------------------------------------------------%
switch opts.Em
    case {'default','quartic'}
        prop.Em = @(l,dp) polyval([3.6567e-14,-2.5074e-10,6.2787e-07,...
            -0.00069181,0.30026],l); % quartic fit to Jellison and Hodgson
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
    case 'Jellison'
        prop.Em_data = getfield(load('Em_Fe_Jellison.mat'),'Em_data');
        prop.Em_gi = griddedInterpolant(prop.Em_data(:,1),prop.Em_data(:,2),'pchip');
        prop.Em = @(l,dp) prop.Em_gi(l);
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
    case 'Hodgson'
        prop.Em_data = getfield(load('Em_Fe_Hodgson.mat'),'Em_data');
        prop.Em_gi = griddedInterpolant(prop.Em_data(:,1),prop.Em_data(:,2),'pchip');
        prop.Em = @(l,dp) prop.Em_gi(l);
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
    case 'Drude' % Jellison and Lowndes? / J. Menser
        prop.omega_p = 2.644E+16;
        prop.tau = 2.385E-16;
        prop.Em = @prop.drude;
        prop.CEmr = 1;
        prop.Emr = @(l1,l2,dp) prop.CEmr.*prop.Em(l1)./prop.Em(l2);
        prop.Eml = @(dp) prop.Em(prop.l_laser,dp);
        
        % Other parameters from Menser, unsure what they are for
        %parameters from Li and Fauchet (1987)
        %omega_p = 2.50e16;
        %tau = 2.12e-16;

        %parameters from Kawamura et al. 
        % omega_p = 2.64E+16;
        % tau = 2.10E-16;
end


prop.opts = opts;

end

