
% Combine the measurement at several wavelengths to reduce noise

function [s,l] = filter_combinewave(s,l,rr)

nn = length(l);
ind0 = 1:rr:nn;
ind1 = rr:rr:nn;

l = mean([l(ind0);l(ind1)]);
for ii=1:length(ind0)
    s(:,:,ii) = mean(s(:,:,ind0(ii):ind1(ii)),3);
end

end

