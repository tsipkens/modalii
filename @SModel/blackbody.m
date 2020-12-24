
% BLACKBODY  Vector implementation of Planck's law.
% AUTHOR: Timothy Sipkens
%=========================================================================%

function Ib = blackbody(T,l)
%-------------------------------------------------------------------------%
% Inputs:
%   l	Wavelength, [nm]
%   T	Temperature, [K]
%
% Outputs:
%   Ib  Blackbody intensity
%-------------------------------------------------------------------------%

p0 = bsxfun(@times,T,reshape(l.*1e-9,1,1,[]));
p1 = 1./(exp(0.0143877696./(p0))); % Number is phi=h*c/kB
    % Wein's approximation used in forward model for scaling factor paper
p2 = 1.191042868e-16./(l.*1e-9).^5; % Number is 2*h*c^2

Ib = bsxfun(@times,p1,reshape(p2,1,1,[]));

end