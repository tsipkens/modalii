function [] = filter_ave_rand(signal,nn)

if nargin<2;
    nn = 20;
end
nn = ceil(length(signal.data(1,:,1))./nn); % modify to number of signals after average

data_ua = signal.data;
signal.data = [];

ind_vec = randperm(length(data_ua(1,:,1)));
for ii=1:nn;
    signal.data(:,ii,:) = nanmean(data_ua(:,ind_vec(ii:nn:end),:),2);
end

end