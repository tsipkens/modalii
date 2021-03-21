
function s = filter_switchneg(s)

bw_zero = 10; % Previously equal to one

s = -s;
for ii=1:length(s(1,1,:)) % Loop through wavelengths
    for jj=1:length(s(1,:,1)) % Loop through signals
        s(:,jj,ii) = s(:,jj,ii)-mean(s(1:bw_zero,jj,ii));
    end
end

end

