
%== SETUP_B ======================================================%
%   Converts b to an appropriate form used elsewhere in the class
%   AUTHOR:  Timothy Sipkens, 2019-07-18

function [bout, Lb, sb] = setup_b(b, variance)

if ~isa(b,'function_handle')
    nn1 = size(b,1);
    nn2 = size(b,2);

    switch variance
        case {'covariance','dependent'}
            Gb = nancov(real(b'))./(nn2-1);
            sb = 1 ./ diag(Gb);
            Gb_i = inv(Gb);
            Lb = chol(stats.Gb_i);
            bout = nanmean(b,2);
        case {'default','variance','independent'}
            sb = nanstd(b,0,2)./sqrt(nn2-1);
            Gb_i = spdiags(1./sb.^2,0,nn1,nn1);
            Lb = chol(Gb_i);
            bout = nanmean(b,2);
        case 'specified'
            bout = nanmean(b,2);
        case 'single-shot'
            sb = (0.05*mean(b)).*...
                ones(nn1,nn2);
            Gb_i = spdiags(1./sb.^2,...
                0,nn1,nn1);
            Lb = sqrt(Gb_i);
            bout = b;
    end

    stats.Gb_i = @(x) Gb_i;
    stats.Lb = @(x) Lb;
    stats.b = @(x) bout;
else
    nn1 = @(x)size(b(x),1);
    nn2 = @(x)size(b(x),2);

    switch variance
        case 'covariance'
            Gb = @(x) nancov(real(b(x)'))./(nn2(x)-1);
            sb = @(x) 1 ./ diag(Gb(x));
            Gb_i = @(x) inv(Gb(x));
            Lb = @(x) chol(Gb_i);
            bout = @(x) nanmean(b(x),2);
        case {'default','variance','independent'}
            sb = @(x) nanstd(b(x),0,2)./sqrt(nn2(x)-1);
            Gb_i = @(x) spdiags(1./sb(x).^2,0,nn1(x),nn1(x));
            Lb = @(x) chol(Gb_i(x));
            bout = @(x) nanmean(b(x),2);
        case 'specified'
            bout = @(x) nanmean(b(x),2);
            Lb = [];
            sb = [];
        case 'single-shot'
            sb = @(x) (0.05*mean(b(x))).*...
                ones(nn1(x),nn2(x));
            Gb_i = @(x) spdiags(1./sb(x).^2,...
                0,nn1(x),nn1(x));
            Lb = @(x) sqrt(Gb_i(x));
            bout = b;
    end
end

end