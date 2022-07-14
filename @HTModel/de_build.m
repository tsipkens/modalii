
% DE_BUILD  Builds function handle for governing equation used in evaluation.
% Currently builds the function using strings and `eval`.
% AUTHOR: Timothy Sipkens, 2018-11-28
% 
% NOTE: The dTdt portion of the equation is initially built up as a string, 
%   with the output depending on htmodel.opts. Allows for quicker
%   evaluation by exlcuding insignificant terms. 
%
% OUTPUTS:
%   dTdt    Function for the change in temperature
%   dmdt    Function for the change in particle mass
%=========================================================================%

function [htmodel, dTdt, dmdt] = de_build(htmodel)

prop = htmodel.prop; % make local copy of material properties
dp = @(mp,T) 1e9.*(6.*mp./(pi.*prop.rho(T))).^(1/3); % output in nm


%-- Mass component of the ODE --------------------------------------------%
dmdt = @(t,T,mp,X) -htmodel.J_evap(prop, T, dp(mp,T));
htmodel.dmdt = dmdt;


%-- Temperature component of the ODE -------------------------------------%
aa = '@(t,T,mp,X)('; % temperature component of the ODE


%-- Conduction model -----------------------------------------------------%
switch htmodel.opts.cond  % default is free molecular regime
    case 'exp'
        htmodel.opts.deMethod = 'none';
    case 'none'
        % Don't add anything.
    otherwise  % free molecular regime
        aa = [aa,'-htmodel.q_cond(prop,T,dp(mp,T))'];
end
%-------------------------------------------------------------------------%


%-- Evaporation model ----------------------------------------------------%
switch htmodel.opts.evap %  Default is free molecular regime
    case 'none'
        % Don't add anything.
    case 'mult'  % free molecular regimes, multiple species
        aa = [aa,'-htmodel.q_evapm(prop,T,dp(mp,T))'];
    otherwise  % free molecular regime, one species
        aa = [aa,'-htmodel.q_evap(prop,T,dp(mp,T))'];
end
%-------------------------------------------------------------------------%


%-- Radiative model ------------------------------------------------------%
switch htmodel.opts.rad % Default is none
    case {'none'}
        % Don't add anything (default).
    otherwise
        aa = [aa,'-htmodel.q_rad(prop,T,dp(mp,T))'];
end
%-------------------------------------------------------------------------%


%-- Abrsorption model ----------------------------------------------------%
switch htmodel.opts.abs % Default is none
    case {'none'}
        % Don't add anything.
    otherwise
        aa = [aa,'+htmodel.q_abs(prop,t,dp(mp,T))'];
end
%-------------------------------------------------------------------------%


%-- Annealing model ------------------------------------------------------%
switch htmodel.opts.ann % Default is none
    case {'include','Michelsen'}
        aa = [aa,'+htmodel.q_ann_Mich(prop,T,dp(mp,T),X)'];
        htmodel.dXdt = @(t,T,mp,X) htmodel.dXdt_fun(@htmodel.q_ann_Mich,prop,T,dp(mp,T),X);
    case {'Sipkens'}
        aa = [aa,'+htmodel.q_ann_Sip(prop,T,dp(mp,T),X)'];
        htmodel.dXdt = @(t,T,mp,X) htmodel.dXdt_fun(@htmodel.q_ann_Sip,prop,T,dp(mp,T),X);
    otherwise
        htmodel.dXdt = @(t,T,mp,X) 0;
end
%-------------------------------------------------------------------------%


%-- Finish and evaluate --------------------------------------------------%
aa = [aa,')./(prop.cp(T).*mp)'];
dTdt = eval(aa);
htmodel.dTdt = dTdt;
%-------------------------------------------------------------------------%


end

