function [W,S_E,S_var,nn,mm,S_var_var] = calcWeight(S)

% scale the data so maximum is 100
S = S./max(mean(S,2)).*100;

% calculate mean / expected value
S_E = mean(S,2);

% calculate variance
S_var = std(S,[],2).^2;

nn = length(S(1,:));
mm = length(S_E);

% calculate standard variance of the variance
S_var_var = ((nn-1).*((nn-1).*moment(S',4)-(nn-3).*(moment(S',2).^2))./(nn^3))';

% modify weight to account for point density
% [ki,xi] = ksdensity(S_E,'bandwidth',1.94);
% % [ki,xi] = ksdensity(S_E);
% Ki = griddedInterpolant(xi,ki);
% Ki_val = Ki(S_E);
% kappa = mean(S_var_var)./mean(S_var_var.*Ki_val); % scale, so mean is unity
% % kappa = exp((-1/mm).*sum(log(Ki_val))); % scale, so product is unity
% S_var_var = S_var_var.*(Ki_val.*kappa);

% weighting matrix or inverse covariance matrix
W = spdiags(1./S_var_var,0,mm,mm);

end

