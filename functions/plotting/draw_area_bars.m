function draw_area_bars(x, individuals, barVals, errVals, coll, faceAlpha, hideZeros)
%DRAW_AREA_BARS Bar, upper error bar and per-brain dots for each area.
%   DRAW_AREA_BARS(x, individuals, barVals, errVals, coll, faceAlpha, hideZeros)
%   draws in the current axes, for area jj at x(jj): a bar of height
%   barVals(jj), an upward error bar of errVals(jj), and one dot per brain
%   (individuals(jj, :)). With hideZeros, zero and NaN values get no dot.

for jj = 1:numel(x)
    ind = individuals(jj, :);
    if hideZeros
        ind(isnan(ind) | ind == 0) = 999;        % outside the y-limits: no dot
    end
    s = bar(x(jj), barVals(jj), 'FaceColor', coll, 'BarWidth', 0.4, ...
        'EdgeColor', 'none', 'FaceAlpha', faceAlpha);
    uistack(s, 'down');
    hold on;
    e = errorbar(x(jj), barVals(jj), errVals(jj), 'CapSize', 4, 'YNegativeDelta', []);
    e.Color = [0.1 0.1 0.1];
    e.LineWidth = 1.0;
    e.XNegativeDelta = [];
    sw = swarmchart(repmat(x(jj), numel(ind), 1), ind', 7, 'filled', ...
        'MarkerFaceColor', coll, 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', coll .* coll .* coll);
    sw.XJitterWidth = 0.1;
    uistack(sw, 'top');
end
end
