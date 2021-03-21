
function s = filter_ave(s,nn)

if nargin<2
    nn = length(s(1,:,1)); % if not specified average all
end
nn = ceil(length(s(1,:,1))./nn); % modify to number of signals after average

s_loc = s;
s = [];

for ii=1:nn
    s(:,ii,:) = nanmean(s_loc(:,ii:nn:end,:),2);
end

end