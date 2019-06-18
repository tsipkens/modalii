function out = blackbody(T,l)
% BLACKBODY Evaluates Planck's law for a given temperature and wavelength.
% l is wavelength in nm
% T is temperature in Kelvin

p0 = bsxfun(@times,T,reshape(l.*1e-9,1,1,[]));
p1 = 1./(exp(0.0143877696./(p0))); % Number is phi=h*c/kB
    % Wein's approximation used in forward model for scaling factor paper
p2 = 1.191042868e-16./(l.*1e-9).^5; % Number is 2*h*c^2
out = bsxfun(@times,p1,reshape(p2,1,1,[]));

end