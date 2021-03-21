
function [s,t,cut] = filter_trimnoise(s,t,t0,opts)


%-- Parse inputs ----------------------------%
tt = 0.04; % Threshold noise level
te = length(t);
ta = 0;
if nargin>1
    if isfield(opts,'t0')
        t0 = opts.t0;
    end
    if isfield(opts,'tt')
        tt = opts.tt;
    end
    if isfield(opts,'te')
        te = opts.te;
    end
    if isfield(opts,'ta')
        ta = opts.ta;
    end
end


%-- Apply filter ------------------------------------------%
if length(s(1,:,1))==1 % Single signal case
    bw = 10;
    end_ind = te; % Length of t vector
    for ii=1:(end_ind-bw)
        sub_sig = s(ii:(ii+bw),:,:);
        mean_sig = squeeze(mean(sub_sig));
        mid_sig = squeeze(sub_sig(floor(bw/2),1,:));
        nsr = max(abs(mid_sig-mean_sig)./mean_sig);
        if nsr>tt
            end_ind = ii;
            break;
        end
    end
    cut = t0:end_ind;
elseif length(s(1,1,:))==1 % Single wavelength/temp. case
    J_std = std(s,0,2)./sqrt(length(s(:,1,1))-1);
    end_ind = min([te,...
        find((J_std(:,:,1)./nanmean(s(:,:,1),2))>tt,1,'first')]);
    cut = t0:end_ind;
else % Multiple wavelength case (Note: only two presently)
    J_std = std(s,0,2)./sqrt(length(s(:,1,1))-1);
    end_ind = min([te,...
        find((J_std(:,:,1)./nanmean(s(:,:,1),2))>tt,1,'first'),...
        find((J_std(:,:,2)./nanmean(s(:,:,2),2))>tt,1,'first')]);
    cut = t0:end_ind;
end
s = s(cut,:,:);
tt = tt(cut)-tt(cut(1))+ta;
t(1) = tt(1)+1e-10;

end


