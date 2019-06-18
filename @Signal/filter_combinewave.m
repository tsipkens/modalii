function [] = filter_combinewave(signal,rr)
% Combine the measurement at several wavelengths to reduce noise

nn = length(signal.l);
ind0 = 1:rr:nn;
ind1 = rr:rr:nn;

signal.l = mean([signal.l(ind0);signal.l(ind1)]);
for ii=1:length(ind0);
    data(:,:,ii) = mean(signal.data(:,:,ind0(ii):ind1(ii)),3);
end

signal.data = data;

end

