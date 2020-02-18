function jcb = jacob_est(stats,mle,min_fun)

delta = 0.005; % to perturb solution from mle

if nargin<3 % just QoIs case
    min_fun = stats.min_fun;
    vec = mle;
end

for ii=1:length(vec)
    vec_h = vec;
    vec_h(ii) = vec_h(ii).*(1+delta);
    f_h = min_fun(vec_h);
    
    vec_l = vec;
    vec_l(ii) = vec_l(ii).*(1-delta);
    f_l = min_fun(vec_l);
    
    jcb(:,ii) = (f_h-f_l)./(2*vec(ii)*delta); % jcb by central difference
end

min_fun(vec); % used to reset parameters in prop to original values

end

