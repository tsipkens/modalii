function [smp,G,R,s] = cred_mcmc(stats,n,mle,sx_pr,theta)

if nargin<4
    sx_pr = 0.04.*mle;
end

if nargin<5
    vec = mle;
    logpdf = @(vec) -(norm(stats.min_fun(vec)./2,2).^2);
    smp = slicesample(vec,n,'logpdf',logpdf,...
        'width',sx_pr); % sampling algorithm
else
    vec = [mle,theta];
    logpdf = @(vec) -(norm(stats.min_fun_theta(...
        vec(1:length(mle)),vec((length(mle)+1):end))./2,2)).^2;
    smp = slicesample(vec,n,'logpdf',logpdf,...
        'width',[sx_pr,stats.st_pr]); % sampling algorithm
end

s = std(smp)';
G = cov(smp);
R = corr(smp);

end

