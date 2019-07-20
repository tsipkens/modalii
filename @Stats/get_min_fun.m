function [min_fun] = get_min_fun(stats)

model = stats.model;

if stats.opts.bFun==0
    switch stats.opts.minimize
        case {'default','levenberg-marquardt','vector'} % minimization considers vector
            min_fun = @(x) ((stats.Lb')*(model(x)-stats.b))./2;

        case {'fminsearch','scalar'} % minimization considers scalar
            min_fun = @(x) norm(((stats.Lb')*(model(x)-stats.b))./2,2);

        case {'vector-xpr'} % function also has priors introduced, vector
            stats.getPriorFun;
            min_fun = @(x) [((stats.Lb')*(model(x)-stats.b))./2;...
                stats.pr_fun(x)];

        otherwise
            disp('ERROR: Invalid minimization technique!');
            return;
    end
    
else
    switch stats.opts.minimize
        case {'default','levenberg-marquardt','vector'} % minimization considers vector
            min_fun = @(x) ((stats.Lb(x)')*(model(x)-stats.b(x)))./2;

        case {'fminsearch','scalar'} % minimization considers scalar
            min_fun = @(x)norm(((stats.Lb(x)')*(model(x)-stats.b(x)))./2,2);

        case {'vector-xpr'} % function also has priors introduced, vector
            stats.get_prior_fun;
            min_fun = @(x) [((stats.Lb(x)')*(model(x)-stats.b(x)))./2;...
                stats.pr_fun(x)];

        otherwise
            disp('ERROR: Invalid minimization technique!');
            return;
    end
end

stats.min_fun = min_fun;

end

