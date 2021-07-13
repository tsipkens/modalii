
function [mle,jcb] = minimize(x0, like, pr, opts)

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

switch opts.minimize
    case {'levenberg-marquardt','default','vector','vector-xpr','vector-tpr'}
        if isfield(opts,'lb')
            options = optimset('Display','off',...
                'Algorithm','trust-region-reflective');
                    % trust-region-reflective required for bounds
            lb = opts.lb;
            ub = opts.ub;
        else
            options = optimset('Display','off','MaxFunEvals',10000,...
                'Algorithm','levenberg-marquardt');
            lb = [];
            ub = [];
        end
        [mle,~,~,~,~,~,jcb] = lsqnonlin(min_fun,x0,lb,ub,options);
        
    case {'trust-region-reflective','trr'}
        if opts.display==1
            options = optimset('MaxFunEvals',10000,...
                'Algorithm','trust-region-reflective');
        else
            options = optimset('MaxFunEvals',10000,'Display','off',...
                'Algorithm','trust-region-reflective');
        end
        [mle,~,~,~,~,~,jcb] = lsqnonlin(min_fun,x0,[],[],options);
        
    case {'fminsearch','scalar'}
        if opts.display==1
            options = optimset('MaxFunEvals',10000,'Display','iter');
        else
            options = optimset('MaxFunEvals',10000);
        end
        mle = fminsearch(min_fun,x0,options);
        jcb = stats.jacob_est(mle);
        
    otherwise
        disp('Invalid minimization technique.');
        return;
end

end

