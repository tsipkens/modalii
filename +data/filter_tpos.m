
function [s,t] = filter_tpos(s,t)

ind_remove = [false;(t(1:(end-1))>t(2:end))];
s(ind_remove,:,:) = [];
t(ind_remove) = [];

end

