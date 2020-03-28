
function [Tp,tp] = get_peak_temp(s,varargin)

[s,t,l] = data.parse_signal(s);

[~,ind_max] = max(s(:,:,1)); % currently based on first wavelength
ind_med = median(ind_max); % use median over shots 
    % (Note: does not incorporate temporal laser jitter)
tp = t(ind_med); % peak time

if isa(varargin{1},'Prop')
    prop = varargin{1};
    prop.l = l; % ensure correct wavelengths are used
    smodel = SModel(prop,{},{},signal);
elseif isa(varargin{1},'SModel')
    smodel = varargin{1};
else
    disp('Invalid input.');
end

[~,Tp] = smodel.IModel(nanmean(s(ind_med,:,:))); % calculate temperature

end

