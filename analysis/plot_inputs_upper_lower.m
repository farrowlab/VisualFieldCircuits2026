%% Brain-wide inputs to one sSC cell type: upper vs lower visual field figures
% Figures for medial (upper visual field) vs lateral (lower visual field)
% sSC injections of one cell type:
%   1. % of inputs per major brain region (dots: brains; line: median),
%      rank-sum test per region, BH-corrected
%   2. % of inputs per major brain region for each brain (stacked bars)
%   3. number of labelled cells per brain
%   4. % of inputs and enrichment per area, grouped by major brain region
%
% Requires the output of load_data_preprocess_inputs_upper_lower.m and
% mysigstar.m (significance stars) on the path.
% Output (P.output): region_ranksum_<celltype>.csv

clearvars; close all;

P = setup_paths();

%% Parameters
celltype = 'ntsr';
sideNames = {'medial', 'lateral'};

% Area selection for figure 4 (each condition in Med or in Lat): passes
% the inclusion criteria of the preprocessing (Threshold), mean enrichment
% (TT column 15) >= minDensity, mean % of inputs outside the sSC (TT
% column 9) >= minPercent; excluded rows removed.
minDensity = 1.6;
minPercent = 1.11;
excludePrefix = 'SCs';                                          % all sSC layers
excludeAreas  = ["MBsen", "MBmot", "MBsta", "ICc", "ICd", "ICe"];

set(0, 'defaultAxesFontSize', 10);

%% Load data
S = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB', 'ParentalAreaN_v2');
BBB = S.BBB;
regionNames = S.ParentalAreaN_v2(:, 1);          % major brain regions (BBB column 6)
nRegions = numel(regionNames);
A = area_indices(BBB);
colors = load(fullfile(P.resources, 'color_schemes.mat'), 'col', 'col_stackbar');
col = colors.col;

sides = {'Med', 'Lat'};                          % Med = upper, Lat = lower visual field
for s = 1:2
    T(s) = load(fullfile(P.processed, [celltype '_' sides{s} '.mat']), 'TT', 'Threshold', 'tableid');
end
nBrains = arrayfun(@(t) size(t.TT{1, 7}, 2), T);

% Rows not counted in the region totals: sSC layers and midbrain summary rows
regionExclude = false(A.nAreas, 1);
regionExclude([A.SC(:); A.MBsen; A.MBmot; A.MBsta]) = true;

%% 1. % of inputs per major brain region
% Sum of the % of inputs (outside sSC) over the areas of each region that
% pass the inclusion criteria
regionPct = cell(1, 2);
for s = 1:2
    [regionPct{s}, present] = region_percent(T(s), nRegions, regionExclude, true);
end
regs = find(present);

figure('Position', [181.8 261.8 775.2 334.4]);
hold on;
gw = min(0.7, 2 / 3);
x = cell(1, 2);
for i = 1:2
    x{i} = (1:nRegions) - gw/2 + (2*i - 1) * gw / 4;
    med = median(regionPct{i}, 2);
    med(~present) = 0;
    h = scatter(x{i}, med', 150, '_', 'LineWidth', 2);
    h.CData = [0.1 0.1 0.1];
    uistack(h, 'top');
    for g = regs
        sw = swarmchart(repmat(x{i}(g), nBrains(i), 1), regionPct{i}(g, :)', 30, ...
            'filled', 'MarkerFaceColor', col(i, :), 'MarkerFaceAlpha', 0.8);
        sw.XJitterWidth = 0.1;
        uistack(sw, 'bottom');
    end
end

% Rank-sum test medial vs lateral per region, BH-corrected
pRegion = arrayfun(@(g) ranksum(regionPct{1}(g, :), regionPct{2}(g, :)), regs);
qRegion = bh_qvalues(pRegion);
for k = 1:numel(regs)
    g = regs(k);
    yl = max([regionPct{1}(g, :), regionPct{2}(g, :)]) + 3;
    mysigstar(gca, mean([x{1}(g) x{2}(g)]), yl, qRegion(k));
end
writetable(table(regionNames(regs), pRegion(:), qRegion(:), 'VariableNames', {'Region', 'ranksum_p', 'ranksum_q'}), ...
    fullfile(P.output, ['region_ranksum_' celltype '.csv']));

for yLine = [25 50 75]
    plot([0, 8.5], [yLine, yLine], '-.', 'Color', [0.5 0.5 0.5]);
end
set(gca, 'XTick', 1:nRegions, 'XTickLabel', regionNames, 'XTickLabelRotation', 0, 'FontSize', 10, 'YTick', 0:25:100);
set(gca, 'box', 'off', 'color', 'none', 'TickDir', 'out', 'Ticklength', [0.001 0.025], 'FontSize', 10);
ylim([0 80]);
ylabel({'Proportion of input', '(%)'}, 'FontSize', 12);
set(gcf, 'renderer', 'Painters');

%% 2. % of inputs per major brain region for each brain
% All areas of each region (no inclusion criteria)
perBrain = [];
for s = 1:2
    [rp, present] = region_percent(T(s), nRegions, regionExclude, false);
    rp(isnan(rp)) = 0;
    perBrain = [perBrain, rp(1:find(present, 1, 'last'), :)]; %#ok<AGROW>
end

figure('Position', [603.67 301.67 1002 362]);
bs = bar(perBrain', 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.6);
hold on;
for i = 1:sum(present)
    bs(1, i).CData = colors.col_stackbar(i, :);
    bs(1, i).FaceAlpha = 0.7;
end
ylim([0 100]);
nAll = sum(nBrains);
set(gca, 'XTick', 1:nAll, 'XTickLabel', num2cell([1:nBrains(1), 1:nBrains(2)]), 'YTick', 0:25:100);
set(gca, 'box', 'off', 'color', 'none', 'TickDir', 'out');
ax = gca;
ax.XAxis.TickLength = [0 0];
ylabel({'% of total inputs'}, 'FontSize', 12);
xlabel('Animal', 'Position', [-0.5, -2.7 -1], 'FontSize', 12);
ax.XRuler.TickLabelGapOffset = +5;
legend(regionNames, 'location', 'northeastoutside', 'AutoUpdate', 'off');
legend('boxoff');
groupRange = {[1, nBrains(1)], [nBrains(1) + 1, nAll]};
for s = 1:2
    r = groupRange{s};
    plot([r(1) - 0.4, r(2) + 0.4], [-12 -12], 'LineWidth', 0.7, 'Clipping', 'off', 'Color', 'k');
    text((r(1) - 0.4 + r(2) + 0.4) / 2, -16, sideNames{s}, 'Clipping', 'off', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
end
set(gca, 'OuterPosition', [0.1 0.1 0.9 0.9]);
set(gcf, 'renderer', 'Painters');

%% 3. Number of labelled cells per brain
% NOTE: sums column 7 over all rows of the area table, so cells in
% overlapping rows (sSC, midbrain summary rows, IC subdivisions) are
% counted more than once.
cellsPerBrain = arrayfun(@(t) sum(cell2mat(t.TT(:, 7)), 1), T, 'UniformOutput', false);
for s = 1:2
    for n = 1:numel(cellsPerBrain{s})
        fprintf('# %s inputs animal %d: %d\n', sides{s}, n, cellsPerBrain{s}(n));
    end
end

figure('Position', [560.2 229 238.4 328]);
xPos = [1 1.5];
for s = 1:2
    h = scatter(xPos(s), mean(cellsPerBrain{s}), 150, '_', 'LineWidth', 2);
    h.CData = [0.1 0.1 0.1];
    uistack(h, 'top');
    hold on;
    for j = 1:numel(cellsPerBrain{s})
        sw = swarmchart(xPos(s), cellsPerBrain{s}(j), 30, 'filled', ...
            'MarkerFaceColor', col(s, :), 'MarkerFaceAlpha', 0.8);
        sw.XJitterWidth = 0.1;
        uistack(sw, 'bottom');
    end
end
set(gca, 'XTick', xPos, 'XTickLabel', sideNames, 'XTickLabelRotation', 0, 'FontSize', 10);
set(gca, 'box', 'off', 'color', 'none', 'TickDir', 'out', 'Ticklength', [0.001 0.025], 'FontSize', 10);
ylabel('Counts per animal', 'FontSize', 12);
set(gcf, 'renderer', 'Painters');

%% 4. % of inputs and enrichment per area, grouped by major brain region
for s = 1:2
    D(s).threshold = T(s).Threshold;
    D(s).percent   = cell2mat(T(s).TT(:, 9));        % % of inputs outside sSC
    D(s).density   = cell2mat(T(s).TT(:, 15));       % enrichment over baseline density
    D(s).count     = cell2mat(T(s).TT(:, 7));
end
acronyms = string(BBB(:, 1));
excluded = startsWith(acronyms, excludePrefix) | ismember(acronyms, excludeAreas);
keepIdx  = select_areas(D, minDensity, minPercent, excluded);
fprintf('Final kept areas: %d of %d\n', numel(keepIdx), numel(acronyms));
for s = 1:2
    fprintf('Total number of included %s %s input cells: %d\n', celltype, sideNames{s}, sum(D(s).count(keepIdx, :), 'all'));
end
for s = 1:2
    fprintf('Average cell counts per %s brain (%s): %d\n', celltype, sideNames{s}, ...
        round(mean(sum(D(s).count(keepIdx, :), 1))));
end

% Order: by median input % (Med + Lat), then by major brain region
labels = BBB(keepIdx, [1 3 6 8]);
[~, plotSeq] = sortrows(median(D(1).percent(keepIdx, :), 2) + median(D(2).percent(keepIdx, :), 2), 'descend');
[~, seq2] = sortrows(labels(plotSeq, :), 3, 'ascend');
plotSeq = keepIdx(plotSeq(seq2));
sortedLabels = BBB(plotSeq, [1 3 6 8]);
region = cell2mat(sortedLabels(:, 3));
lr = numel(plotSeq);

pctS = cell(1, 2); denS = cell(1, 2);
for s = 1:2
    pctS{s} = D(s).percent(plotSeq, :);
    denS{s} = D(s).density(plotSeq, :);
end
x = cell(1, 2);
for k = 1:2
    x{k} = (1:2:lr*2) - 0.5 + (2*k - 1) / 4;
end

figure('Position', [617 261 765.6 352]);
spacing = 5;
tiledlayout(2*spacing + 1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- % of inputs (bars: median; error bars: MAD; dots: brains) ---
ax1 = nexttile([spacing 1]);
for i = 1:2
    draw_area_bars(x{i}, pctS{i}, median(pctS{i}, 2), mad(pctS{i}, 1, 2), col(i, :), 0.8, false);
end
shade_brain_regions(region, 30);

% % --- Significance per area, medial vs lateral ---
% % uncorrected t-test:
% % pvals = arrayfun(@(jj) ttest2_pvalue(pctS{1}(jj, :), pctS{2}(jj, :)), 1:lr);
% % rank-sum test, Benjamini-Hochberg corrected:
% pvals = arrayfun(@(jj) ranksum(pctS{1}(jj, :), pctS{2}(jj, :)), 1:lr);
% padj = bh_qvalues(pvals);
% for jj = 1:lr
%     xvector = mean([x{1}(jj) x{2}(jj)]);
%     yl = max([pctS{1}(jj, :), pctS{2}(jj, :)]) + 2;
%     mysigstar(gca, xvector, yl, padj(jj));      % sigstars on top
% end

label_brain_regions(region, 26);
ylim([0 30]);
ax = gca;
set(ax, 'XTick', (1:2:lr*2), 'TickLabelInterpreter', 'none', ...
    'XTickLabel', [], 'XTickLabelRotation', 90, 'FontSize', 10, 'FontName', 'Arial', ...
    'GridColor', [.5 .5 .5], 'GridLineStyle', '--', 'XGrid', 'off', 'color', 'none', 'TickDir', 'out', ...
    'YGrid', 'on', 'YTick', 0:10:100, 'YTickLabel', string(0:10:100), ...
    'GridAlpha', 0.5, 'GridLineWidth', 0.8);
ax.XAxis.TickLength = [0 0.5];
ylabel('% of input');
curX = xlim;
xlim([0 curX(2)]);
box off;
legend('', 'Location', 'northeast');
legendNames = ["up", "low"];
for j = 1:2
    patch([NaN NaN NaN NaN], [NaN NaN NaN NaN], col(j, :), 'FaceAlpha', 1, ...
        'EdgeColor', 'none', 'DisplayName', legendNames(j));
end
hold off;

% --- Area names ---
ax2 = nexttile;
text((1:2:lr*2), repmat(0.5, 1, lr), sortedLabels(:, 1), 'HorizontalAlignment', 'center', ...
    'FontSize', 9, 'Rotation', 90);
set(gca, 'color', 'none');
axis off;

% --- Enrichment over baseline density (bars: median; error bars: MAD) ---
ax3 = nexttile([spacing 1]);
for i = 1:2
    draw_area_bars(x{i}, denS{i}, median(denS{i}, 2), mad(denS{i}, 1, 2), col(i, :), 0.8, false);
end
yAuto = ylim;
shade_brain_regions(region, yAuto(2));

% % --- Significance per area, medial vs lateral ---
% % uncorrected t-test:
% % pvals = arrayfun(@(jj) ttest2_pvalue(denS{1}(jj, :), denS{2}(jj, :)), 1:lr);
% % rank-sum test, Benjamini-Hochberg corrected:
% pvals = arrayfun(@(jj) ranksum(denS{1}(jj, :), denS{2}(jj, :)), 1:lr);
% padj = bh_qvalues(pvals);
% for jj = 1:lr
%     xvector = mean([x{1}(jj) x{2}(jj)]);
%     yl = max([denS{1}(jj, :), denS{2}(jj, :)]) + 2;
%     mysigstar(gca, xvector, yl, padj(jj));      % sigstars on top
% end

ylim([0 300]);
ax = gca;
set(ax, 'box', 'off', 'XTick', (1:2:lr*2), 'TickLabelInterpreter', 'none', ...
    'XTickLabel', [], 'YTickLabelRotation', 0, 'FontSize', 10, 'FontName', 'Arial', ...
    'GridColor', [.5 .5 .5], 'GridLineStyle', '--', 'XGrid', 'off', 'color', 'none', ...
    'YDir', 'reverse', 'TickDir', 'out', ...
    'YGrid', 'on', 'YTick', 0:100:350, 'YTickLabel', string(0:100:350), ...
    'GridAlpha', 0.5, 'GridLineWidth', 0.8);
ax.XAxis.TickLength = [0 0.5];
ax.XAxis.Visible = 'off';
ylabel('relative density (D/mean D)');
curX = xlim;
xlim([0 curX(2)]);
box off;

linkaxes([ax1 ax2 ax3], 'x');
set(gcf, 'renderer', 'Painters');
set(gca, 'color', 'none');


%% Local functions
function [perRegion, present] = region_percent(T, nRegions, excludeRows, applyThreshold)
% Sum of the % of inputs (TT column 9) over the areas of each major brain
% region (TT column 4), per brain. With applyThreshold, areas that fail the
% inclusion criteria count as 0. Regions without areas are NaN.
pct = cell2mat(T.TT(:, 9));
if applyThreshold
    pct(T.Threshold == 0, :) = 0;
end
region = cell2mat(T.TT(:, 4));
present = ismember(1:nRegions, region);
perRegion = nan(nRegions, size(pct, 2));
for g = find(present)
    perRegion(g, :) = sum(pct(region == g & ~excludeRows, :), 1);
end
end

function p = ttest2_pvalue(a, b)
% p-value of a two-sample t-test (for the commented significance options)
[~, p] = ttest2(a, b);
end
