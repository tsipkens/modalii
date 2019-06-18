function [mle,jcb] = minimize(stats,x0,opts)

if or(nargin<2,isempty(x0))
    x0 = ones(length(htmodel.x),1);
end

if nargin<3;
    opts = stats.opts;
end

min_fun = stats.getMinFun;

% Add check to make sure that model output and data are the same size

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
        if opts.display==1;
            options = optimset('MaxFunEvals',10000,...
                'Algorithm','trust-region-reflective');
        else
            options = optimset('MaxFunEvals',10000,'Display','off',...
                'Algorithm','trust-region-reflective');
        end
        [mle,~,~,~,~,~,jcb] = lsqnonlin(min_fun,x0,[],[],options);
        
    case {'fminsearch'}
        if opts.display==1;
            options = optimset('MaxFunEvals',10000,'Display','iter');
        else
            options = optimset('MaxFunEvals',10000);
        end
        mle = fminsearch(min_fun,x0,options);
        jcb = stats.jcbEst(mle);
        % [mle,~,~,~,jcb] = fminunc(min_fun,x0,options);
        
    case 'GN'
        L=chol(inv(diag(stats.sb(:)).^2));
        z=L*data(:);
        Opts.uiprint=0;
        [mle{1},mle{2}]=GNiteration(@GNimp,z,x0,Opts,model,L);
        jcb = [];
        
    otherwise
        disp('Invalid minimization technique.');
        return;
end

end

