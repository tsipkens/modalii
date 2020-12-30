
% DE_GEN  Builds function handle for governing equation used in evaluation.
% Currently builds the function using strings and `eval`.
% AUTHOR: Timothy Sipkens, 2018-11-28
% 
% NOTE: The dTdt portion of the equation is initially built up as a string, 
%   with the output depending on htmodel.opts. Allows for quicker
%   evaluation by exlcuding insignificant terms. 
%
% INPUTS:
%   dp0     Vector of nanoparticle diameters for which ODEs are to be solved, [nm]
%
% OUTPUTS:
%   Tout    Time-resolved temperature with a single column per input diameter, [K]
%   dpo     Time-resolved nanoparticle diameter, same format as above, [nm]
%   mpo     Time-resolved nanoparticle mass, same format as above, [fraction]
%   Xo      Time-resolved anneealed fraction, same format as above, [fraction]
%=========================================================================%

function [htmodel, dTdt, dmdt] = de_gen(htmodel)

prop = htmodel.prop; % make local copy of material properties
dp = @(mp,T) 1e9.*(6.*mp./(pi.*prop.rho(T))).^(1/3); % output in nm


%-- Mass component of the ODE --------------------------------------------%
dmdt = @(t,T,mp,X) -htmodel.J_evap(T,dp(mp,T));
htmodel.dmdt = dmdt;

tools.textheader('Generating heat transfer ODE');


%-- Temperature component of the ODE -------------------------------------%
aa = '@(t,T,mp,X)('; % temperature component of the ODE


%-- Conduction model -----------------------------------------------------%
switch htmodel.opts.cond % Default is free molecular regime
    case 'exp'
        disp('  CONDUCTION:  Exponential');
        htmodel.opts.deMethod = 'none';
    case 'none'
        disp('  CONDUCTION:  -');
        % Do nothing.
    otherwise % Free molecular regime
        disp('  CONDUCTION:  Free molecular');
        aa = [aa,'-htmodel.q_cond(T,dp(mp,T))'];
end
%-------------------------------------------------------------------------%


%-- Evaporation model ----------------------------------------------------%
switch htmodel.opts.evap % Deafult is free molecular regime
    case 'none'
        disp('  EVAPORATION: -');
        % Do nothing.
    case 'mult' % Free molecular regimes, multiple species
        disp('  EVAPORATION: Multiple species');
        aa = [aa,'-htmodel.q_evapm(T,dp(mp,T))'];
    otherwise % Free molecular regime, one species
        disp('  EVAPORATION: Free molecular');
        aa = [aa,'-htmodel.q_evap(T,dp(mp,T))'];
end
%-------------------------------------------------------------------------%


%-- Radiative model ------------------------------------------------------%
switch htmodel.opts.rad % Default is none
    case {'include'}
        disp('  RADIATION:   Incl.');
        aa = [aa,'-htmodel.q_rad(T,dp(mp,T))'];
    otherwise
        % Do nothing.
end
%-------------------------------------------------------------------------%


%-- Abrsorption model ----------------------------------------------------%
switch htmodel.opts.abs % Default is none
    case {'include'}
        disp('  ABSORPTION:  Incl.');
        aa = [aa,'+htmodel.q_abs(t,dp(mp,T))'];
    otherwise
        disp('  ABSORPTION:  -');
        % Do nothing. 
end
%-------------------------------------------------------------------------%


%-- Annealing model ------------------------------------------------------%
switch htmodel.opts.ann % Default is none
    case {'include','Michelsen'}
        disp('  ANNEALING:   Michelsen');
        aa = [aa,'+htmodel.q_ann_Mich(T,dp(mp,T),X)'];
        htmodel.dXdt = @(t,T,mp,X) htmodel.dXdt_fun(@htmodel.Qann_Mich,T,dp(mp,T),X);
    case {'Sipkens'}
        disp('  ANNEALING:   Sipkens');
        aa = [aa,'+htmodel.q_ann_Sip(T,dp(mp,T),X)'];
        htmodel.dXdt = @(t,T,mp,X) htmodel.dXdt_fun(@htmodel.Qann_Sip,T,dp(mp,T),X);
    otherwise
        htmodel.dXdt = @(t,T,mp,X) 0;
end
%-------------------------------------------------------------------------%


%-- Finish and evaluate --------------------------------------------------%
aa = [aa,')./(prop.cp(T).*mp)'];
dTdt = eval(aa);
htmodel.dTdt = dTdt;
%-------------------------------------------------------------------------%

tools.textheader();

end


