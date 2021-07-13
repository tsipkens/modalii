
% JACOB_EST  Function to estimate the Jacobian using perturbation analysis. 
%  
%  JCB = stats.jacob_est(MLE,LIKE) uses the likelihood function LIKE and
%  maximum likelihood estaimte MLE, without prior information and using the
%  default options.
%  
%  JCB = stats.jacob_est(MLE,LIKE,PR) adds a prior function specified by PR
%  (expects a maximum a posteriori estimate for MLE). 
%  
%  JCB = stats.jacob_est(MLE,LIKE,PR,OPTS) adds an options structure.
%  
%  AUTHOR: Timothy Sipkens, 2021-04-22

function jcb = jacob_est(mle, like, pr, opts)

%-- Parse inputs ----------------------------------------------%
if ~exist('pr', 'var'); pr = []; end
if isempty(pr); pr = @(x) []; end

if ~exist('opts', 'var'); opts = []; end
if isempty(opts); opts = struct(); end
if ~isfield(opts, 'minimize'); opts.minimize = 'default'; end

if or(strcmp(opts.minimize, 'fminsearch'), ...
        strcmp(opts.minimize, 'scalar'))
    min_fun = @(x) norm([like(x); pr(x)]);
else
    min_fun = @(x) [like(x); pr(x)];
end
%--------------------------------------------------------------%

delta = 0.005; % to perturb solution from mle

for ii=1:length(mle)
    mle_h = mle;
    mle_h(ii) = mle_h(ii).*(1+delta);
    f_h = min_fun(mle_h);
    
    mle_l = mle;
    mle_l(ii) = mle_l(ii).*(1-delta);
    f_l = min_fun(mle_l);
    
    jcb(:,ii) = (f_h-f_l)./(2*mle(ii)*delta); % jcb by central difference
end

min_fun(mle); % used to reset parameters in prop to original values

end

