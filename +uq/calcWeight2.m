function [W,S_E,S_var,nn,mm,S_var_var] = calcWeight2(S)

% scale the data so maximum is 100
S = S./max(mean(S,2)).*100;

% calculate mean / expected value
S_E = mean(S,2);

% calculate variance
S_var = std(S,[],2).^2;

nn = length(S(1,:));
mm = length(S_E);

% calculate variance of the variance by sectioning
% oo = 40;
% S_var_set = zeros(mm,oo);
% for ii=1:oo;
%     S_var_set(:,ii) = std(S(:,ii:oo:end),[],2).^2;
% end
S_var_var = (S_var.^2).*(2/(nn-1));

% weighting matrix or inverse covariance matrix
W = spdiags(1./S_var_var,0,mm,mm);

end

