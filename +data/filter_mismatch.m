
% will also determine whether peak is unreasonably late

function s = filter_mismatch(s,tp)

tp_m = ceil(mean(mean(tp)));

% Find max index and check if great than tp+n for any channel
% If true, remove shot. 
ii_elim = [];
for ii=2:length(s(1,:,1))
    if or(max2(s(:,ii,1))>(tp_m+10),...
            max2(s(:,ii,2))>(tp_m+10))
        ii_elim = [ii_elim,ii];
    end
end
s(:,ii_elim,:)=[];
tp(:,ii_elim) = [];

disp(['Removed ',num2str(length(ii_elim)),' signals due to peak mismatch.']);

% Find shots zero voltage at peak
% If true, remove shot.
ii_zero = [];
for ii=2:length(s(1,:,1))
    if or((s(tp_m,ii,1))==0,...
            (s(tp_m,ii,2))==0)
        ii_zero = [ii_zero,ii];
    end
end
s(:,ii_zero,:)=[];
tp(:,ii_zero) = [];

disp(['Removed ',num2str(length(ii_zero)),' signals due to zero entries.']);

tp = tp_m; % Update and store tp
t0 = tp_m;

end




function ind = max2(s,dim)

if nargin>1
    [~,ind] = max(s,dim);
else
    [~,ind] = max(s);
end

end

