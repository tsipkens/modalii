
% FMODEL    Evaluates forward model (T > J) using given E(m) and temperature decay.
% Author:   Timothy Sipkens
%=========================================================================%

function [Jout] = FModel(smodel,T,Em)

p3 = Em(smodel.l,smodel.prop.dp0)./(smodel.l.*1e-9)./smodel.data_sc;

J = smodel.blackbody(T,smodel.l); % evaluate Planck's Law

Jout = bsxfun(@times,J,reshape(p3,1,1,[]));
    % multiply Planck's Law by necessary values to return incandescence

end

