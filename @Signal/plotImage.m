function [] = plotImage(signal,ind)

if nargin<2
    ind = 1;
end

imagesc(signal.l,signal.t,squeeze(signal.data(:,ind,:)));

end

