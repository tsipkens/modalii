function [Tp,tp,F0] = getPeakTemp(signal,varargin)

F0 = signal.F0;

[~,ind_max] = max(signal.data(:,:,1)); % currently based on first wavelength
ind_med = median(ind_max); % use median over shots 
    % (Note: does not incorporate temporal laser jitter)
tp = signal.t(ind_med); % peak time

if isa(varargin{1},'Prop')
    prop = varargin{1};
    prop.l = signal.l; % ensure correct wavelengths are used
    smodel = SModel(prop,{},{},signal);
elseif isa(varargin{1},'SModel')
    smodel = varargin{1};
else
    disp('Invalid input.');
end

[~,Tp] = smodel.IModel(nanmean(signal.data(ind_med,:,:))); % calculate temperature

end

