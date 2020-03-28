
function s = filter_ave_rand(s,nn)

if nargin<2
    nn = 20; % number of effective signals to output
end
nn = ceil(length(s(1,:,1))./nn);
    % modify to number of signals after average

s_ua = s; % unaveraged signals
s = [];

ind_vec = randperm(length(s_ua(1,:,1)));
for ii=1:nn
    s(:,ii,:) = nanmean(s_ua(:,ind_vec(ii:nn:end),:),2);
end

end