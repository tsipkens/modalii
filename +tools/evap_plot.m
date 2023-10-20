
% EVAP_PLOT  Plots components of one or more evaporation models. 
%  This includes vapor pressure, latent heat, and mass of vapor.
%  
%  AUTHOR: Timothy Sipkens, 2023-10-19

function [] = evap_plot(prop)

if ~iscell(prop); prop = {}; end

Tvec = 1e3:1:5e3;

for ii=1:length(prop)
    subplot(3, 1, 1);
    plot(Tvec, prop{ii}.hv(Tvec, prop{ii}));
    hold on;
    ylabel('Latent heat, hv');
    
    subplot(3, 1, 2);
    plot(Tvec, prop{ii}.pv(Tvec, 100, [], prop{ii}));  % evaluate at 100 nm (~ avoids size effects)
    hold on;
    ylabel('Vapor pressure, pv');

    subplot(3, 1, 3);
    if isa(prop{ii}.Mv, 'function_handle')
        plot(Tvec, prop{ii}.Mv(Tvec, prop{ii}) .* 1e3);
    else
        plot(Tvec, prop{ii}.Mv .* ones(size(Tvec)) .* 1e3);
    end
    hold on;
    ylabel('Vapor molar mass, Mv');
end

for jj=1:3
    subplot(3, 1, jj);
    hold off;

    if jj==2
        yline(101325);
        set(gca, 'YScale', 'log');
    end

    if jj < 3
        h = gca; h.XAxis.TickLabels = {};
    end
end

end
