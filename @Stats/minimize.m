function [mle,jcb] = minimize(stats,x0,opts)

if or(nargin<2,isempty(x0))
    x0 = ones(length(htmodel.x),1);
end

if nargin<3
    opts = stats.opts;
end

min_fun = stats.get_min_fun; % generate function to minimize

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
        
    case {'fminsearch'}
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

