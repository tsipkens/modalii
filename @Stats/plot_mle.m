function [] = plot_mle(stats,mle,t)

if nargin < 3
    t = 1:length(stats.data(mle));
end

clf;
if stats.opts.bFun==0 % handle if b is a function
    b0 = stats.b;
else
    b0 = stats.b(mle);
end
t_pl = t-t(1);

if 1 % plot data error bounds
    hold on;
    % plot(t-t(1),[b0+2.*stats.sb(mle),b0-2.*stats.sb(mle)]);   
    fill([t_pl;flipud(t_pl)],...
        [b0-2.*stats.sb(mle);flipud(b0+2.*stats.sb(mle))],...
        [0.94 0.94 0.94],'LineStyle','none');
%     area([t_pl;flipud(t_pl)],...
%         [b0-2.*stats.sb(mle);flipud(b0+2.*stats.sb(mle))]);  
    hold off;
end
hold on;
plot(t_pl,b0,'.','markers',4,'Color',[0 0.443 0.737]);
plot(t_pl,stats.model(mle),'k');
hold off;

end

