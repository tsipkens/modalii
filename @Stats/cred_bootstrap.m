function [smp,G,R,s] = cred_bootstrap(stats,n,mle)

b_str = stats.b;
res0 = (stats.model(mle)-stats.b(mle))./stats.sb(mle);

smp = zeros(n,length(mle));
for ii=1:n
    res_ii = datasample(res0,length(res0));
    T_ii = stats.model(mle)+(res_ii.*stats.sb(mle));
    stats.b = @(varargin) T_ii;
    smp(ii,:) = stats.minimize(mle);
    if mod(ii,25)==0
        disp([num2str(round(ii/n*100)),'% completed.']);
    end
end

s = std(smp)';
G = cov(smp);
R = corr(smp);

stats.b = b_str;

end

