
% Q_ABS  Rate of laser energy input from the laser. 
% 
%  INPUTS:
%   t   Time, [ns]
%   dp	Nanoparticle diameter, [nm]
%
%  OUTPUTS:
%   q	Rate of laser energy uptake by the nanoparticle, [W]
%  
%  AUTHOR: Timothy Sipkens, 2018-12-17

function [q,Cabs] = q_abs(htmodel, prop, t, dp)

dp = dp .* 1e-9; % convert to meters so everything is in SI units
tlp = prop.tlp * 1e-9; % convert from ns to s
tlm = prop.tlm * 1e-9; % convert from ns to s
t = t .* 1e-9; % convert from ns to s

F1 = (prop.F0') .* (100*100); % laser fluence, [J/m2]
    % Note: the above includes conversion from [J/cm2]

% Evaluate absorption cross section.
Cabs = pi ^ 2 .* (dp .^ 3) ./ (prop.l_laser*1e-9) .* ...
    prop.Eml(dp * 1e9);


%-- Setup laser profile as function handle, f(t) -------------------------%
switch htmodel.opts.abs
    case {'tophat','square'} % square laser profile
        f = @(t) F1.*(heaviside(t - (tlm - tlp/2))...
            -heaviside(t - (tlm + tlp/2))) ./ (tlp);
        
    case {'Gaussian','normal','include'} % Gaussian temporal laser profile
        sigma = tlp ./ (2 * sqrt(2 * log(2))); % makes tlp the FWHM of the pulse
        
        % tlm = sigma; % --> changed for model selection work?
        % was sigma = 5e-9
        
        f = @(t) F1 .* normpdf(t, tlm, sigma);
            
    case {'lognormal'}
        Sk = 0.9282; % based on logn fit to Michelsen, 2007 profile
        f_A = nthroot(Sk^2+sqrt(Sk^4+4*Sk^2)+2,3)/nthroot(2,3);
        f_s = sqrt(log(f_A+1/f_A-1));
        f_m = log(tlp / (2 * sqrt(2 * log(2))) / ...
            sqrt(exp(f_s ^ 2) - 1)) - 1/2 * f_s ^ 2;
        f = @(t) F1 .* pdf('Lognormal', t + exp(f_m), f_m, f_s);
            % corresponds to lognormal with E(X)=3*tlp, Var(X) = 8*(tlp^2)*ln^2(2)
            
end
%-------------------------------------------------------------------------%

q = Cabs .* f(t); % evaluate rate of energy uptake

end

