
% BLACKBODY  Vector implementation of Planck's law.
%  Static function, i.e., evaluation does not rely on or use SModel
%  properties.
%  
%  IB = SModel.blackbody(T, L) computes blackbody radiation for a
%  temperature of T in K and wavelength of L in nm. The temperature can be
%  a 2D matrix. Output, IB, will have the same size as T in the first two 
%  dimensions and information across the different wavelengths in the third
%  dimension. 
%  
%   AUTHOR: Timothy Sipkens

function Ib = blackbody(T,l)

p0 = bsxfun(@times, T, reshape(l .* 1e-9, 1, 1, []));

% Wein's approximation used in forward model for scaling factor paper
p1 = 1 ./ (exp(0.0143877696 ./ (p0))); % Number is phi=h*c/kB
p2 = 1.191042868e-16 ./ (l .* 1e-9) .^ 5; % constant givne by 2*h*c^2

Ib = bsxfun(@times, p1, reshape(p2, 1, 1, []));

end