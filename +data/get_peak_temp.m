
% GET_PEAK_TEMP  Estimate peak temperature from experimental signals.
% AUTHOR: Timothy Sipkens, 2020-01-02
%=========================================================================%

function [Tp,tp] = get_peak_temp(s, varargin)

%-- Parse inputs ---------------------------------------------------------%
if isa(varargin{1}, 'SModel')
    smodel = varargin{1};
    prop = smodel.prop;
elseif isa(varargin{1}, 'struct')
    prop = varargin{1};
    prop.l = l; % ensure correct wavelengths are used
    smodel = SModel(prop, {}, {}, signal);
else
    disp('Invalid input.');
end


[s,t,l] = data.parse_signal(s);


[~,ind_max] = max(s(:,:,1)); % currently based on first wavelength
ind_med = median(ind_max); % use median over shots 
    % (Note: does not incorporate temporal laser jitter)
tp = t(ind_med); % peak time

[~, Tp] = smodel.IModel(prop, nanmean(s(ind_med,:,:))); % calculate temperature

end

