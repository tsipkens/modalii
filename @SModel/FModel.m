function [Jout] = FModel(smodel,T,Em)
% FMODEL Evaluates incnadescence using given E(m) and temperature decay.

p3 = Em(smodel.l,smodel.prop.dp0)./(smodel.l.*1e-9)./smodel.data_sc;
J = smodel.blackbody(T,smodel.l);
Jout = bsxfun(@times,J,reshape(p3,1,1,[]));

end

