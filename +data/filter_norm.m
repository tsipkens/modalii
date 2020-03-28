
function s = filter_norm(s)

for ii=1:length(s(1,:,1)) % Loop through shots
    maxsig = max(max(s(:,ii,:)));
    s(:,ii,:) = s(:,ii,:)./maxsig;
end

end

