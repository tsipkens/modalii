
% CALCRATIOPYROMETRY    Evaluate temperature by two-colour pyrometry.
% Author:               Timothy Sipkens
%=========================================================================%

function [Tout,Cout,s_T,out] = calcRatioPyrometry(smodel,J1,J2)

l = smodel.l;
Emr = smodel.prop.Emr(l(1),l(2),smodel.prop.dp0); % two-colour pyrometry
Jr = J1./J2;


%-- Basic ratio calculation ----------------------------------------------%
Tout = (0.0143877696*(1/(l(2)*1e-9)-1/(l(1)*1e-9)))./...
    log(Jr.*(((l(1)/l(2))^6)/Emr)); % phi=0.0143...
Tout = real(Tout);
    % Note: Imaginary values result from negative values in J matrices
s_T = std(Tout,[],2);


%-- Calculate scaling constant -------------------------------------------%
Cout = bsxfun(@rdivide,J2,...
    (smodel.blackbody(Tout,l(end)).*smodel.prop.Em(l(end),smodel.prop.dp0)./...
    (l(end)*1e-9.*smodel.data_sc)));
out = [];
out.s_C = std(Cout,[],2);


%-- Calculate correlation ------------------------------------------------%
ntime = length(Tout(:,1));
r_TC = zeros(ntime,1);
for ii=1:ntime
    r_TC(ii) = corr(Cout(ii,:)',Tout(ii,:)');
end
out.r_TC = r_TC;

end

