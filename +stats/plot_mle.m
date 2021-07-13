
% PLOT_MLE  Plot the MLE solution and data.
%  
%  NOTE: Smoothing is applied to data intervals.

function [] = plot_mle(mle, A, b, t, sb)

if ~isa(b, 'function_handle'); b = @(x) b; end

if ~exist('t', 'var'); t = []; end
if isempty(t); t = 1:length(b(mle)); end

if ~exist('sb', 'var'); sb = []; end
if isempty(sb); sb = 0 .* b(mle); end
if isa(sb, 'function_handle'); sb = sb(mle); end

clf;
b0 = b(mle);
t_pl = t(:) -t(1);

if 1 % plot data error bounds
    hold on;
    sb_mod = smooth(sb, 0.07);
    fill([t_pl; flipud(t_pl)], ...
        [b0 - 2.*sb_mod; flipud(b0 + 2.*sb_mod)], ...
        [0.9 0.9 0.9], 'LineStyle', 'none');
    % area([t_pl;flipud(t_pl)],...
    %   [b0-2.*stats.sb(mle);flipud(b0+2.*stats.sb(mle))]);  
    hold off;
end

hold on;
plot(t_pl, b0, '.', 'markers', 5, ...
    'Color', [0 0.443 0.737]);
plot(t_pl, A(mle), 'k');
hold off;

end

