function [] = filter_ave(signal,nn)

if nargin<2;
    nn = length(signal.data(1,:,1)); % if not specified average all
end
nn = ceil(length(signal.data(1,:,1))./nn); % modify to number of signals after average

data_loc = signal.data;
signal.data = [];

for ii=1:nn;
    signal.data(:,ii,:) = nanmean(data_loc(:,ii:nn:end,:),2);
end

signal.opts.isave = true;

end