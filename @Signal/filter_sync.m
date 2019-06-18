function [] = filter_sync(signal,varargin)
% automatically aligns peaks of different signals

sz = size(signal.data);
ind_shift = zeros(1,sz(2),sz(3));
sync_channels = false; % sync channels shot by shot

% Parse inputs *******************
ii=1;
while ii<=(nargin-1);
    if strcmp(varargin{ii},'shift')
        ind_shift(1,1,:) = varargin{ii+1};
        ii = ii+2; % Skip an input
    elseif strcmp(varargin{ii},'sync_channels')
        sync_channels = varargin{ii+1};
        ii = ii+2; % Skip an input
    else
        disp(['Input ''',varargin{ii},''' is not an option.']);
        ii = ii+2; % Skip an input
    end
end

% Determine maximum for shift **************
ind_end = length(signal.data(:,1,:));
[~,ind_max] = max(signal.data); % Find index of max in time
ind_sync = ind_max;
if sync_channels==0; % Only use first channel
    ind_sync = repmat(ind_sync(:,:,1),[1,1,length(ind_sync(1,1,:))]);
end

ind_sync = ind_sync-min(min(ind_sync))+1;
ind_sync = ind_sync+ind_shift; % Incorporate added shift
ind_absmax = max(max(max(ind_sync))); % Maximum shifted entry
temp = zeros(ind_end-ind_absmax,length(signal.data(1,:,1)),...
    length(signal.data(1,1,:)));
for jj=1:length(temp(1,:,1)); % shot loop
    for ii=1:length(temp(1,1,:)); % wavelength loop
        temp(:,jj,ii) = signal.data(ind_sync(1,jj,ii):(ind_end-ind_absmax+ind_sync(1,jj,ii)-1),jj,ii);
    end
end
signal.data = temp;

signal.t = signal.t(1:(ind_end-ind_absmax));

end

