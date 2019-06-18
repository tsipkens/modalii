function [bout,Lb_i,sb,Gb_i] = getNoiseInfo(stats,bin)

if stats.opts.bFun==0 % handle case that b is a function
    nn1 = size(bin,1);
    nn2 = size(bin,2);
    
    switch stats.opts.variance
        case {'covariance','dependent'}
            Gb_i = inv(nancov(real(bin'))./...
                (nn2-1));
            Lb_i = chol(Gb_i);
            bout = nanmean(bin,2);
        case {'default','variance','independent'}
            sb = nanstd(bin,0,2)./sqrt(nn2-1);
            Gb_i = spdiags(1./sb.^2,0,nn1,nn1);
            Lb_i = chol(Gb_i);
            bout = nanmean(bin,2);
        case 'specified'
            bout = nanmean(bin,2);
        case 'single-shot'
            sb = (0.05*mean(bin)).*...
                ones(nn1,nn2);
            Gb_i = spdiags(1./sb.^2,...
                0,nn1,nn1);
            Lb_i = sqrt(Gb_i);
            bout = bin;
    end
else
    nn1 = @(x)size(bin(x),1);
    nn2 = @(x)size(bin(x),2);
    
    switch stats.opts.variance
        case 'covariance'
            Gb_i = @(x) inv(nancov(real(bin(x)'))./...
                (nn2(x)-1));
            Lb_i = @(x) chol(Gb_i);
            bout = @(x)nanmean(bin(x),2);
        case {'default','variance','independent'}
            sb = @(x) nanstd(bin(x),0,2)./sqrt(nn2(x)-1);
            Gb_i = @(x) spdiags(1./sb(x).^2,0,nn1(x),nn1(x));
            Lb_i = @(x) chol(Gb_i(x));
            bout = @(x) nanmean(bin(x),2);
        case 'specified'
            bout = @(x) nanmean(bin(x),2);
            Lb_i = [];
            sb = [];
            Gb_i = [];
        case 'single-shot'
            sb = @(x) (0.05*mean(bin(x))).*...
                ones(nn1(x),nn2(x));
            Gb_i = @(x) spdiags(1./sb(x).^2,...
                0,nn1(x),nn1(x));
            Lb_i = @(x) sqrt(Gb_i(x));
            bout = bin;
    end
end

end

