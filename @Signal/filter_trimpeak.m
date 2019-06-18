function [] = filter_trimpeak(signal)

[~,ind_max] = max(signal.data(:,:,1)); % currently based on first wavelength
ind_med = median(ind_max); % use median over shots 
    % (Note: does not incorporate temporal laser jitter)
signal.t = signal.t(ind_med:end);
signal.data = signal.data(ind_med:end,:,:);

end

