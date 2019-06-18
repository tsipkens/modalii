function [] = filter_trimnoise(signal,opts)

% Handle options **********************************************************
t0 = evalt0(signal); % Starting trim time
tt = 0.04; % Threshold noise level
te = length(signal.t);
ta = 0;
if nargin>1;
    if isfield(opts,'t0')
        t0 = opts.t0;
    end
    if isfield(opts,'tt')
        tt = opts.tt;
    end
    if isfield(opts,'te');
        te = opts.te;
    end
    if isfield(opts,'ta');
        ta = opts.ta;
    end
end
signal.t0 = t0; % Update the initial value

% Apply filter ************************************************************
if length(signal.data(1,:,1))==1 % Single signal case
    bw = 10;
    end_ind = te; % Length of t vector
    for ii=1:(end_ind-bw);
        sub_sig = signal.data(ii:(ii+bw),:,:);
        mean_sig = squeeze(mean(sub_sig));
        mid_sig = squeeze(sub_sig(floor(bw/2),1,:));
        nsr = max(abs(mid_sig-mean_sig)./mean_sig);
        if nsr>tt;
            end_ind = ii;
            break;
        end
    end
    signal.cut = signal.t0:end_ind;
elseif length(signal.data(1,1,:))==1 % Single wavelength/temp. case
    J_std = std(signal.data,0,2)./sqrt(length(signal.data(:,1,1))-1);
    end_ind = min([te,...
        find((J_std(:,:,1)./nanmean(signal.data(:,:,1),2))>tt,1,'first')]);
    signal.cut = signal.t0:end_ind;
else % Multiple wavelength case (Note: only two presently)
    J_std = std(signal.data,0,2)./sqrt(length(signal.data(:,1,1))-1);
    end_ind = min([te,...
        find((J_std(:,:,1)./nanmean(signal.data(:,:,1),2))>tt,1,'first'),...
        find((J_std(:,:,2)./nanmean(signal.data(:,:,2),2))>tt,1,'first')]);
    signal.cut = signal.t0:end_ind;
end
signal.data = signal.data(signal.cut,:,:);
signal.t = signal.t(signal.cut)-signal.t(signal.cut(1))+ta;
signal.t(1) = signal.t(1)+1e-10;

end

function [t0] = evalt0(signal)

if isempty(signal.t0)
    t0 = 15; % origianlly 20
else
    t0 = signal.t0; % Use default stored in signal
end

end

