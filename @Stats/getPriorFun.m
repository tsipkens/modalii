function [pr_fun,pr_fun_theta] = getPriorFun(stats)

switch stats.opts.prior
    case {'default','none','uninformed'} % do not use a prior
        pr_fun = @(x) 1;
        
    case {'independent'} % only apply priors on QoIs
        pr_fun = @(x) ((x-stats.x_pr)./stats.sx_pr)';
        
    case {'xpr','normal'} % only apply priors on QoIs
        pr_fun = @(x) (stats.Lx_ipr')*(x-stats.x_pr)';
        
    otherwise
        disp('ERROR: Invalid prior scheme!');
        return;
end

stats.pr_fun = pr_fun;

end

