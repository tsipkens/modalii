function [] = filter_tpos(signal)

ind_remove = [false;(signal.t(1:(end-1))>signal.t(2:end))];
signal.data(ind_remove,:,:) = [];
signal.t(ind_remove) = [];

end

