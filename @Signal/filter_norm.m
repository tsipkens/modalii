function [] = filter_norm(signal)

for ii=1:length(signal.data(1,:,1)) % Loop through shots
    maxsig = max(max(signal.data(:,ii,:)));
    signal.data(:,ii,:) = signal.data(:,ii,:)./maxsig;
end

end

