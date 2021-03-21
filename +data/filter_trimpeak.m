
function [s,t] = filter_trimpeak(s,t)

[~,ind_max] = max(s(:,:,1)); % currently based on first wavelength
ind_med = median(ind_max); % use median over shots 
    % (Note: does not incorporate temporal laser jitter)
t = t(ind_med:end);
s = s(ind_med:end,:,:);

end

