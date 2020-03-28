
% OUTLIER  Remove outliers from two color LII data using Thompson tau.
% Examines each channel in sequence - if one element is removed, the other 
% is also removed
% Author: Timothy Sipkens
%=========================================================================%

function s = outlier(s)

ntot = numel(s);
T = real(s);
 
T(T < 0) = NaN(1);

alph = 0.975;
tau = @(N)icdf('T',alph,N-2).*(N-1)./(sqrt(N).*sqrt(N-2+icdf('T',alph,N-2).^2));
    % from t-distribution
    
done = false;
ll = 0;
while (done == false)
    done = true;
    mu = nanmean(T,2);
    sigma = nanstd(T,[],2);
    n = sum(not(isnan(T)),2);
    n0 = size(T);
  
    % Determine the value of S times Tau
    TS=tau(n).*sigma; % for samples of size 3 < n < 40
   
    delta = abs(T - repmat(mu,[1,n0(2),1]));
    [deltamax,ii] = max(delta,[],2);
    if any(any(deltamax > TS))
        kk = repmat((1:n0(1))',[1,1,n0(3:end)]);
        ii(deltamax <= TS)=[]; % determine sample no. to remove
        kk(deltamax <= TS)=[]; % determine time to remove
        
        jj = false(size(T));
        jj(sub2ind(size(jj),kk,ii))=true;
        jj = repmat(any(jj,3),[1,1,n0(3:end)]);
        
        T(jj) = NaN(1);
        done = false;
        ll = ll+sum(sum(jj));
    end
end

% signal.data;
% signal.Gamma_inv = inv(nancov(T'));
% signal.data = mu;
% signal.sigma = sigma./sqrt(n-1);

disp(['Removed ',num2str(ll(1)),' of ',num2str(ntot)...
    ' points, as they are outliers.']);

end
