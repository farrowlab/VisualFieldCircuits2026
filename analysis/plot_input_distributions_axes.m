%% Distribution of input cells along the brain axes: upper vs lower visual field
% For one cell type, kernel density of the position of the labelled input
% cells along the M-L, A-P and D-V axes, for medial (upper visual field,
% top row) and lateral (lower visual field, bottom row) sSC injections.
% Lines: individual brains; grey area: all cells of all brains pooled.
% Cells in the sSC and in areas failing the inclusion criteria are left out.
% Medial and lateral distributions are compared by Cohen's d, the
% Kolmogorov-Smirnov D statistic and their overlap (all on pooled cells).
%
% Requires the output of load_data_preprocess_inputs_upper_lower.m.
% Output (P.output): DensityPlot_<celltype>_<date>.eps,
%                    distribution_stats_<celltype>.csv

clearvars; close all;

P = setup_paths();

%% Parameters
celltype   = 'grp';
axesNames  = {'ML', 'AP', 'DV'};
axisColumn = struct('AP', 1, 'DV', 2, 'ML', 3);   % columns of the coordinates (AP, DV, ML)
smoothSpan = 5;                                    % moving-average span of the density curves
yLimits    = [0 0.7];

set(0, 'defaultAxesFontSize', 10);

%% Load cell coordinates per brain
colors = load(fullfile(P.resources, 'color_schemes.mat'), 'col');
col = colors.col;

sides = {'Med', 'Lat'};                            % Med = upper, Lat = lower visual field
brainXYZ = cell(1, 2);
for s = 1:2
    T = load(fullfile(P.processed, [celltype '_' sides{s} '.mat']), 'TT', 'Threshold', 'CoordMerge');
    keep = ~startsWith(string(T.TT(:, 1)), 'SCs') & T.Threshold ~= 0;
    nBrains = size(T.CoordMerge, 2);
    brainXYZ{s} = cell(1, nBrains);
    for i = 1:nBrains
        brainXYZ{s}{i} = cell2mat(T.CoordMerge(keep, i));    % cells x [AP DV ML]
    end
end

%% Plot densities and compare medial vs lateral
figure('WindowState', 'maximized', 'Position', [488.2 255.4 651.2 406.4]);
stats = table('Size', [numel(axesNames), 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, ...
    'VariableNames', {'Axis', 'CohensD', 'KS_D', 'Overlap'});

for a = 1:numel(axesNames)
    axName = axesNames{a};
    c = axisColumn.(axName);
    pooled = cell(1, 2);
    sub = gobjects(1, 2);

    for s = 1:2
        sub(s) = subplot(2, 3, a + 3 * (s - 1));
        hold on;

        % Individual brains
        pooled{s} = [];
        for i = 1:numel(brainXYZ{s})
            v = brainXYZ{s}{i}(:, c);
            [f, xj] = ksdensity(v, 'Support', kde_support(v, axName));
            h = plot(xj, smooth(f, smoothSpan)', 'Color', col(s, :));
            h.LineWidth = 0.7;
            pooled{s} = [pooled{s}; v];
        end

        % All cells pooled
        [Y, xi] = ksdensity(pooled{s}, 'Support', kde_support(pooled{s}, axName));
        Y = smooth(Y, smoothSpan)';
        ar = area(xi(Y >= 0), Y(Y >= 0));
        ar.FaceColor = [0.8 0.8 0.8];
        ar.EdgeColor = 'none';
        uistack(ar, 'bottom');

        set(gca, 'color', 'none', 'box', 'off', 'FontSize', 11, 'TickDir', 'out');
    end

    % Axis ticks and labels
    allTicks = [sub(1).XTick sub(2).XTick];
    minVal = floor(min([pooled{1}; pooled{2}]));
    maxVal = ceil(max([pooled{1}; pooled{2}]));
    switch axName
        case 'AP'      % AP was made positive in the preprocessing: label as -mm from bregma
            for s = 1:2
                sub(s).XTick = -2:2:7;
                sub(s).XTickLabel = split(num2str(sub(s).XTick * (-1)));
                sub(s).XLim = [minVal maxVal];
                sub(s).XTickLabelRotation = 0;
            end
            sub(2).XLabel.String = 'A-P (Bregma: mm)';
        case 'DV'
            for s = 1:2
                sub(s).XTick = 1:2:7;
                sub(s).XTickLabel = split(num2str(sub(s).XTick * (-1)));
                sub(s).XTickLabelRotation = 0;
            end
            sub(1).XLim = [minVal maxVal];
            sub(2).XLabel.String = 'D-V (Surface: mm)';
        otherwise      % ML
            for s = 1:2
                sub(s).XTick = min(allTicks) + 1:2:max(allTicks);
                sub(s).XTickLabelRotation = 0;
            end
            sub(1).XLim = [minVal maxVal];
            sub(2).XLabel.String = 'M-L (Middleline: mm)';
    end
    sub(2).XLabel.FontSize = 12;
    sub(1).YLim = yLimits;
    sub(2).YLim = yLimits;
    linkaxes(sub, 'xy');

    if a == 1
        lab = ylabel(sub(2), 'Probability density', 'FontSize', 14);
        pos = get(lab, 'Position');
        lab.Position(1) = pos(1) - 0.55;
        lab.Position(2) = pos(2) * 2.5;
    end

    % Medial vs lateral
    [d, ksD, overlap] = distribution_difference(pooled{1}, pooled{2});
    stats(a, :) = {axName, d, ksD, overlap};
    fprintf('\n%s axis\n', axName);
    fprintf('  Cohen''s d = %.3f\n', d);
    fprintf('  KS D statistic = %.3f  (%s)\n', ksD, ks_size_label(ksD));
    fprintf('  Distribution overlap = %.3f (1 = identical, 0 = no overlap)\n', overlap);
end

set(gcf, 'renderer', 'Painters');
print(gcf, fullfile(P.output, sprintf('DensityPlot_%s_%s.eps', celltype, ...
    datestr(datetime('today'), 'yyyymmdd'))), '-depsc', '-r300');
writetable(stats, fullfile(P.output, ['distribution_stats_' celltype '.csv']));


%% Local functions
function sup = kde_support(v, axName)
% Bounded support for ksdensity: D-V starts at the brain surface (0)
if strcmp(axName, 'DV')
    sup = [0 ceil(max(v) + 0.01)];
else
    sup = [floor(min(v) - 0.01) ceil(max(v) + 0.01)];
end
end

function [d, ksD, overlap] = distribution_difference(x1, x2)
% Cohen's d (x1 - x2), two-sample KS statistic, and overlap of the two
% kernel densities (shared area / union area)
pooledStd = sqrt((std(x1)^2 + std(x2)^2) / 2);
d = (mean(x1) - mean(x2)) / pooledStd;
[~, ~, ksD] = kstest2(x1, x2);

[f1, xi1] = ksdensity(x1);
[f2, xi2] = ksdensity(x2);
xi = linspace(min([xi1, xi2]), max([xi1, xi2]), 1000);
f1 = interp1(xi1, f1, xi, 'linear', 0);
f2 = interp1(xi2, f2, xi, 'linear', 0);
overlap = sum(min(f1, f2)) / sum(max(f1, f2));
end

function s = ks_size_label(D)
if D < 0.1
    s = 'negligible';
elseif D < 0.3
    s = 'small';
elseif D < 0.5
    s = 'medium';
else
    s = 'large';
end
end
