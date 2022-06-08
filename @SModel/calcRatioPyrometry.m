
% CALCRATIOPYROMETRY  Evaluate temperature by two-colour pyrometry.
%  
%  T = SModel.calcRatioPyrometry(J1, J2) evaluates the temperature using
%  ratio pyrometry and the optical properties embedded in SModel.
%  
%  T = SModel.calcRatioPyrometry(J1, J2, EMR) adds an input for the ratio
%  of the absorption function at the two wavelengths. If excluded, the
%  value is extracted from SModel.prop. 
%  
%  T = SModel.calcRatioPyrometry(J1, J2, EMR, IDX) allows for selecting
%  specific wavelength indices on which to peform pyrometry (e.g., first 
%  and second wavelength, IDX = [1,2]).
%  
%  AUTHOR: Timothy Sipkens, 2017

function [To, Co, s_T, out] = calcRatioPyrometry(smodel, J1, J2, Emr, idx)

%-- Parse inputs ---------------------------------------------------------%
l = smodel.l;  % local copy

if length(l) > 2  % check the number of wavelengths
    if ~exist('idx', 'var'); idx = []; end
    if length(idx) ~= 2
        error(['More than two wavelengths in SModel. ', ...
            'Ratio pyrometry not performed.']);
    end
    l = l(idx);
end

if ~exist('Emr', 'var'); Emr = []; end
if isempty(Emr)
    Emr = smodel.prop.Emr(l(1), l(2), smodel.prop.dp0); % ratio of Em at two wavelengths
end
%-------------------------------------------------------------------------%

Jr = J1 ./ J2;  % ratio of incandescence

%-- Basic ratio calculation ----------------------------------------------%
To = (0.0143877696 * (1/(l(2)*1e-9) - 1/(l(1)*1e-9))) ./ ...
    log(Jr .* (((l(1)/l(2))^6) / Emr)); % phi=0.0143...

% Note: Imaginary values result from negative values in J matrices. 
% This avoids them. 
To = real(To);

% Compute the range of temperatures if multiple signals were given.
s_T = std(To,[],2);


%-- Calculate scaling constant -------------------------------------------%
%   NOTE: Does not account for annealing.
Co = bsxfun(@rdivide,J2,...
    (smodel.blackbody(To,l(end)) .* smodel.prop.Em(l(end),smodel.prop.dp0,0) ./ ...
    (l(end)*1e-9.*smodel.data_sc)));
out = [];
out.s_C = std(Co,[],2);


%-- Calculate correlation ------------------------------------------------%
ntime = length(To(:,1));
r_TC = zeros(ntime,1);
for ii=1:ntime
    r_TC(ii) = corr(Co(ii,:)',To(ii,:)');
end
out.r_TC = r_TC;

end

