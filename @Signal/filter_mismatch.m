function [] = filter_mismatch(signal,tp)
% will also determine whether peak is unreasonably late

if nargin>1;
    signal.tp = tp;
elseif isempty(signal.tp)
    tp = zeros(length(signal.data(1,1,:)),length(signal.data(1,:,1)));
    for jj=1:length(tp(:,1));
        for ii=1:length(tp(1,:));
            tp(jj,ii) = max2(signal.data(:,ii,jj));
        end
    end
else
    tp = signal.tp;
end
tp_m = ceil(mean(mean(tp)));

% Find max index and check if great than tp+n for any channel
% If true, remove shot. 
ii_elim = [];
for ii=2:length(signal.data(1,:,1));
    if or(max2(signal.data(:,ii,1))>(tp_m+10),...
            max2(signal.data(:,ii,2))>(tp_m+10));
        ii_elim = [ii_elim,ii];
    end
end
signal.data(:,ii_elim,:)=[];
tp(:,ii_elim) = [];

disp(['Removed ',num2str(length(ii_elim)),' signals due to peak mismatch.']);

% Find shots zero voltage at peak
% If true, remove shot.
ii_zero = [];
for ii=2:length(signal.data(1,:,1));
    if or((signal.data(tp_m,ii,1))==0,...
            (signal.data(tp_m,ii,2))==0);
        ii_zero = [ii_zero,ii];
    end
end
signal.data(:,ii_zero,:)=[];
tp(:,ii_zero) = [];

disp(['Removed ',num2str(length(ii_zero)),' signals due to zero entries.']);

signal.tp = tp_m; % Update and store tp
signal.t0 = tp_m;

end

function ind = max2(data,dim)

if nargin>1;
    [~,ind] = max(data,dim);
else
    [~,ind] = max(data);
end

end

