
function prop = exper_apb17_Mo(prop,~)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = props.init(); end

%-- Conduction properties ------------------------------------------------%
prop.Tg = 298;
prop.Pg = 101325;
prop.Ti = prop.Tg;
prop.tlp = 8;
prop.tlm = 30;
prop.l_laser = 1064;

%-- Particle size and signal properties ----------------------------------%
prop.dp0 = 20;
prop.sigma = 0; % Default monodisperse

end

