
function prop = x_cenide_ge(prop, ~)

if ~exist('prop', 'var'); prop = []; end
if isempty(prop); prop = props.init(); end

%-- Conduction properties ------------------------------------------------%
prop.Tg = 1560;
prop.Pg = 1e4;

%-- Particle size and signal properties ----------------------------------%
prop.dp0 = 20;
prop.sigma = 0; % Default monodisperse

end
