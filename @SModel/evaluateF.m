
% EVALUATEF  Evaluate spectroscopic forward model (T/htmodel -> J). 
%  
%  JOUT = SModel.evaluateF(X) evaluates the spectroscopic model for
%  a given set of QoI, X. Temperature if computed via the embedded heat 
%  transfer model. Output is a scaled incandescence, JOUT.
%  
%  AUTHOR: Timothy Sipkens

function [Jo, mp] = evaluateF(smodel, x)

htmodel = smodel.htmodel; % embedded heat transfer model

if and(isa(smodel.prop,'struct'),~isfield(smodel.prop,'C_J'))
    smodel.prop.C_J = 1;
end

%-- Update x values in prop struct ---------------------------------------%
[smodel, prop] = tools.update_prop(smodel, x);
[htmodel, ~] = tools.update_prop(htmodel, x);
%-------------------------------------------------------------------------%


%-- POLYDISPERSE ---------------------------------------------------------%
%   In this case, integrate incandescence over the size distribution.
if prop.sigma > 0.005 % currently models with lognormal
    dp0 = prop.dp0; % store current dp0 and apply as geometric mean
    
    % Discrete particle size in logspace to evaluate distribution.
    dp_1 = logninv(0.01, real(log(prop.dp0)), prop.sigma);
    dp_2 = logninv(0.9999, real(log(prop.dp0)), prop.sigma);
    dp_x_diff = (dp_2-dp_1)/50;  % element sizes
    dp_x = (dp_1:dp_x_diff:dp_2)';  % differentiated space
    
    % Compute the weight. 
    w = (lognpdf(dp_x, log(dp0), prop.sigma) .* dp_x_diff .* ...
        pi .* (dp_x.*1e-9).^3)';
    w = w ./ sum(w);  % normalize the weight
    
    % Evaluate temperature and incandescence
    [T, ~, mp, X] = htmodel.de_solve(prop, dp_x);  % solve differential equation
    mpr = real(mp ./ mp(1,:));
    
    Jo = bsxfun(@times, w, ...
        smodel.FModel(prop, T, prop.Em, X)); % evaluate forward model for J
    Jo = Jo .* mpr;
    Jo = sum(Jo, 2);
    Jo = Jo .* prop.C_J; % scale incandescence for stability
    mp = []; % not currently calculating mass ratio for distribution
    

%-- MONODISPERSE ---------------------------------------------------------%
%   If the distribution is narrow enough, skip integration over 
%   size distribution and evaluate temperature directly (much faster).
else % for monodisperse case, simply evaluate the ODE directly
    [T, ~, mp, X] = htmodel.de_solve(prop, prop.dp0); % solve heat transfer model at dp0
    
    Jo = smodel.FModel(prop, T, prop.Em, X); % evaluate forward model for J
    Jo = Jo .* prop.C_J; % scale incandescence by corresponding factor
    
    % Incorporate mass loss (i.e., a reduction in volume fraction).
    mpr = real(mp ./ mp(1));
    Jo = bsxfun(@times, Jo, mpr);
    
end
%-------------------------------------------------------------------------%

end
