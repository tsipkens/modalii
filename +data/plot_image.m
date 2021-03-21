
function [] = plot_image(s,t,l,ind)

if nargin<2
    ind = 1;
end

imagesc(l,t,squeeze(s(:,ind,:)));

end

