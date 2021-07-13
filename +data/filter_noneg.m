
% FILTER_NONNEG  Replace negative signals with NaNs.

function s = filter_noneg(s)

s(s < 0) = NaN;

end
