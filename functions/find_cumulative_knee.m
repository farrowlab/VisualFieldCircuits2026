function [optimalThresh, distances, kneeIdx] = find_cumulative_knee(percents, threshRange)
%FIND_CUMULATIVE_KNEE Input-percentage threshold at the knee of the cumulative curve.
%   [optimalThresh, distances, kneeIdx] = FIND_CUMULATIVE_KNEE(percents, threshRange)
%
%   percents     areas-by-brains input percentages
%   threshRange  candidate thresholds (vector)
%
%   For each candidate threshold the % of areas whose mean input % is at or
%   above it is computed. The knee (Kneedle method) is the point of this
%   curve, after normalising both axes to [0, 1], that lies farthest from
%   the straight line joining its first and last points.

pm = mean(percents, 2);
n  = numel(pm);
threshRange = threshRange(:)';

cumulative = sum(pm >= threshRange, 1)' / n * 100;

x = (threshRange' - min(threshRange)) / (max(threshRange) - min(threshRange));
y = (cumulative - min(cumulative)) / (max(cumulative) - min(cumulative));

distances = abs((y(end) - y(1)) * x - (x(end) - x(1)) * y + x(end) * y(1) - y(end) * x(1)) ...
    / sqrt((y(end) - y(1))^2 + (x(end) - x(1))^2);

[~, kneeIdx] = max(distances);
optimalThresh = threshRange(kneeIdx);

fprintf('Search boundary: %.2f%% | Knee threshold: %.2f%% (%d areas retained)\n', ...
    threshRange(end), optimalThresh, sum(pm >= optimalThresh));
end
