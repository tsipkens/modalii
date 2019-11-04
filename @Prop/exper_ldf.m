
% Experimental parameters for a standard laminar diffusion flame (Gulder burner) experiment.

function prop = exper_ldf(prop,~)

%-- Conduction properties --------------------------%
prop.Tg = 1730;
prop.Pg = 101325;
prop.Ti = prop.Tg;
prop.tlp = 8;
prop.tlm = 0;
prop.l_laser = 1064;

%-- Particle size and related properties ----------%
prop.dp0 = 30;
prop.sigma = 0; % Default monodisperse

end

