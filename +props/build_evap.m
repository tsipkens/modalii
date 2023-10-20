
% BUILD_EVAP  A function to combine evaporation of multiple species. 
%  
%  AUTHOR: Timothy Sipkens, 2023-10-12

function [prop] = build_evap(prop)

% Clausius-Clapyeron for partial pressures of each species
for ii=1:length(prop.evap)
    prop.evap(ii).pv = eval(['@(T,dp,hv,prop) 101325.*exp(prop.evap(', num2str(ii), ...
        ').hv(prop.evap(', num2str(ii), ').Tref)./8.3145.*(1./prop.evap(', ...
        num2str(ii), ').Tref-1./T))']);
end

% Compile combined quantities. 
pv_str = '@(T,dp,hv,prop) prop.evap(1).pv(T,[],[],prop)';
for ii=2:length(prop.evap)
    pv_str = [pv_str, ' + prop.evap(', num2str(ii), ').pv(T,[],[],prop)'];
end
prop.pv = eval([pv_str, ';']);

for ii=1:length(prop.evap)
    prop.evap(ii).pv_frac = @(T,dp,hv,prop) prop.evap(ii).pv(T,dp,hv,prop) ./ prop.pv(T,dp,hv,prop);
end

% Combine remaining quantities. 
Mv_str = '@(T,prop) 0';
hv_str = '@(T,prop) (0';
alpham_str = '@(T,prop) (0';
for ii=1:length(prop.evap)
    Mv_str = [Mv_str, ' + prop.evap(', num2str(ii), ').pv_frac(T,[],[],prop)', ...
         ' .* prop.evap(', num2str(ii), ').Mv'];
    hv_str = [hv_str, ' + prop.evap(', num2str(ii), ').pv_frac(T,[],[],prop)', ...
         ' .* prop.evap(', num2str(ii), ').hv(T,prop)'];
    alpham_str = [alpham_str, ' + 0.5.*prop.evap(',num2str(ii) , ...
        ').pv(T,[],[],prop).*prop.evap(', num2str(ii), ').hv(T,prop)./((prop.evap(', num2str(ii), ...
        ').Mv)).^0.5'];
end

prop.Mv = eval([Mv_str, ';']);
prop.mv = @(T,prop) prop.Mv(T,prop) .* 1.660538782e-24;
prop.hv = eval([hv_str, ') ./ prop.Mv(T,prop);']);

alpham_str = [alpham_str, ')./(prop.pv(T,[],[],prop).*prop.Mv(T,prop).^0.5.*prop.hv(T,prop));'];
prop.alpham = eval(alpham_str);

end
