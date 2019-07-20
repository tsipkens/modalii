function [post,x,y] = eval_post(stats,X,Y)

grid = [];
grid(:,:,1) = X;
grid(:,:,2) = Y;

post = zeros(length(grid(:,1,1)),length(grid(1,:,1)));
for ii=1:length(grid(:,1,1))
    for jj=1:length(grid(1,:,1))
        post(ii,jj) = -(norm(stats.min_fun(squeeze(grid(ii,jj,:))),2)^2)/2;
    end
end

end

