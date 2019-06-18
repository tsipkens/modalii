function [] = filter_switchneg(signal)

bw_zero = 10; % Previously equal to one

signal.data = -signal.data;
for ii=1:length(signal.data(1,1,:)) % Loop through wavelengths
    for jj=1:length(signal.data(1,:,1)) % Loop through signals
        signal.data(:,jj,ii) = signal.data(:,jj,ii)-mean(signal.data(1:bw_zero,jj,ii));
    end
end

end

