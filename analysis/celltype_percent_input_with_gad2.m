%% Inputs to Ntsr1+ (wide-field), Grp+ (narrow-field) and Gad2+ sSC neurons, all injections pooled
% Uses the pooled data of each cell type (upper and lower visual field
% injections combined) to
%   1. select the areas to analyse,
%   2. show the input-percentage threshold (knee of the cumulative distribution),
%   3. plot % of inputs per area and cell type, grouped by major brain region,
%      with a Kruskal-Wallis test across cell types per area,
%   4. correlate input % between two cell types,
%   5. test per area for an input bias between ntsr and grp (analyze_input_bias),
%   6. cluster areas by their co-variation across brains, per cell type,
%      and compare the two cluster structures (Mantel test).
%
% Requires the output of load_data_preprocess_inputs_all_celltypes.m.
% Output (P.output): celltype_abbr.csv, kruskalwallis_celltypes.csv,
%                    bias_ntsr_vs_grp_pooled.mat/.csv

clearvars; close all;

P = setup_paths();

%% Parameters
celltypes = {'ntsr', 'grp', 'gad2'};
cellNames = {'WF', 'NF', 'Gad2'};                % names in figure legends and labels

% Area selection, combined across cell types (each condition in at least
% one cell type): passes the inclusion criteria of the preprocessing
% (Threshold), mean enrichment (TT column 15) >= minDensity, mean % of
% inputs outside the sSC (TT column 9) >= minPercent; excluded rows removed.
minDensity = 1.6;
minPercent = 1.11;
excludePrefix = 'SCs';                                          % all sSC layers
excludeAreas  = ["MBsen", "MBmot", "MBsta", "ICc", "ICd", "ICe", "LSd"];

% Knee analysis: cell types used, and the percent criterion applied before
% the knee is computed (minPercent reproduces the current figure)
kneeCelltypes  = [1 2];                         % ntsr, grp
kneeMinPercent = minPercent;
threshRange    = 0:0.01:20;

% Correlation of input % between two cell types (indices into celltypes)
corrPair = [1 2];

% Hierarchical clustering: cut the dendrogram at this fraction of its height
clusterCutoff = 0.62;

nBoot = 10000;            % bootstrap of the correlation
nPerm = 10000;            % Mantel test permutations
seed  = 42;

set(0, 'defaultAxesFontSize', 10);

%% Load data
S   = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB');
BBB = S.BBB;
colors = load(fullfile(P.resources, 'color_schemes.mat'), 'col', 'col_bias_med', 'col_bias_lat');
ctColors = colors.col([5 6 3], :);              % ntsr, grp, gad2

nCt = numel(celltypes);
for c = 1:nCt
    S  = load(fullfile(P.processed, ['OD_' celltypes{c} '_pooled.mat']), 'OD_cell');
    T  = load(fullfile(P.processed, [celltypes{c} '_pooled.mat']), 'TT', 'Threshold', 'tableid');
    D(c).count     = cell2mat(T.TT(:, 7));        % areas x brains
    D(c).percent   = cell2mat(T.TT(:, 9));        % % of inputs outside sSC
    D(c).density   = cell2mat(T.TT(:, 15));       % enrichment over baseline density
    D(c).threshold = T.Threshold;
    D(c).nTotal    = [S.OD_cell.nTotalminusSC{T.tableid}];
end

%% Select areas (same set for all cell types)
acronyms = string(BBB(:, 1));
excluded = startsWith(acronyms, excludePrefix) | ismember(acronyms, excludeAreas);
keepIdx  = select_areas(D, minDensity, minPercent, excluded);
fprintf('Final kept areas: %d of %d\n', numel(keepIdx), numel(acronyms));

labels = BBB(keepIdx, [1 3 6 8]);               % acronym, name, parent group, volume
areaNames = labels(:, 1)';
pct = cell(1, nCt);
cnt = cell(1, nCt);
for c = 1:nCt
    pct{c} = D(c).percent(keepIdx, :);
    cnt{c} = D(c).count(keepIdx, :);
    fprintf('Total number of included %s input cells: %d\n', celltypes{c}, sum(cnt{c}, 'all'));
    fprintf('Average cell counts per %s brain: %d\n', celltypes{c}, round(mean(sum(cnt{c}, 1))));
end

%% Input-percentage threshold: knee of the cumulative distribution
kneeIdx = select_areas(D, minDensity, kneeMinPercent, excluded);
knee = struct();
for k = 1:numel(kneeCelltypes)
    c = kneeCelltypes(k);
    fprintf('\nAnalyzing cell type: %s\n', celltypes{c});
    knee(k).pct = D(c).percent(kneeIdx, :);
    [knee(k).thresh, knee(k).distances, knee(k).idx] = find_cumulative_knee(knee(k).pct, threshRange);
end
kneeMean = mean([knee.thresh]);
fprintf('\nFINAL THRESHOLD: %.2f%%\n', kneeMean);
for k = 1:numel(kneeCelltypes)
    fprintf('  %s: %d areas retained\n', celltypes{kneeCelltypes(k)}, ...
        sum(mean(knee(k).pct, 2) >= kneeMean));
end
plot_knee(knee, kneeMean, threshRange, ctColors(kneeCelltypes, :));

%% % of inputs per area and cell type, grouped by major brain region
% Bars: median across brains; error bars: median absolute deviation; dots: brains
medPct = cell2mat(cellfun(@(x) median(x, 2), pct, 'UniformOutput', false));
madPct = cell2mat(cellfun(@(x) mad(x, 1, 2), pct, 'UniformOutput', false));

% Order: by ntsr median input %, then by major brain region
[~, plotSeq] = sortrows(medPct(:, 1), 'descend', 'MissingPlacement', 'last');
[~, seq2] = sortrows(labels(plotSeq, :), 3, 'ascend');
plotSeq = plotSeq(seq2);
sortedLabels = labels(plotSeq, :);
writetable(cell2table(sortedLabels(:, 1:2)), fullfile(P.output, 'celltype_abbr.csv'));

% Kruskal-Wallis test across cell types per area, BH-corrected
kwP = nan(numel(plotSeq), 1);
for jj = 1:numel(plotSeq)
    vals  = cellfun(@(x) x(plotSeq(jj), :), pct, 'UniformOutput', false);
    group = cellfun(@(v, k) k * ones(1, numel(v)), vals, num2cell(1:nCt), 'UniformOutput', false);
    kwP(jj) = kruskalwallis([vals{:}]', [group{:}]', 'off');
end
kwQ = bh_qvalues(kwP);
writetable(table(sortedLabels(:, 1), kwP, kwQ, 'VariableNames', {'Area', 'KW_p', 'KW_q'}), ...
    fullfile(P.output, 'kruskalwallis_celltypes.csv'));

plot_grouped_percent(cellfun(@(x) x(plotSeq, :), pct, 'UniformOutput', false), ...
    medPct(plotSeq, :), madPct(plotSeq, :), sortedLabels, ctColors, cellNames);

%% Correlation of input % between two cell types
% Per area: median and median absolute deviation of log input % across
% brains; bootstrap over brains (separately per cell type) for the CIs
p1 = pct{corrPair(1)};
p2 = pct{corrPair(2)};
nonzeroVals = [p1(p1 > 0); p2(p2 > 0)];
pseudocount = min(nonzeroVals) / 10;              % below the smallest nonzero value
log1 = log(p1 + pseudocount);
log2v = log(p2 + pseudocount);
med1 = median(log1, 2); med2 = median(log2v, 2);
mad1 = mad(log1, 1, 2); mad2 = mad(log2v, 1, 2);

n1 = size(log1, 2);
n2 = size(log2v, 2);
x_fit  = linspace(min(med1(isfinite(med1))), max(med1(isfinite(med1))), 100);
boot_r = zeros(nBoot, 1);
y_boot = zeros(nBoot, numel(x_fit));
rng(seed);
for b = 1:nBoot
    bv1 = median(log1(:, randi(n1, n1, 1)), 2);
    bv2 = median(log2v(:, randi(n2, n2, 1)), 2);
    valid = isfinite(bv1) & isfinite(bv2);
    if sum(valid) > 2
        boot_r(b) = corr(bv1(valid), bv2(valid), 'type', 'Pearson');
        y_boot(b, :) = polyval(polyfit(bv1(valid), bv2(valid), 1), x_fit);
    else
        boot_r(b) = NaN;
        y_boot(b, :) = NaN;
    end
end
ci_r = prctile(boot_r, [2.5 97.5]);
valid = isfinite(med1) & isfinite(med2);
[r, p] = corr(med1(valid), med2(valid), 'type', 'Pearson');

figure('Position', [351.4 335.4 412 324]);
hold on;
fill([x_fit, fliplr(x_fit)], [prctile(y_boot, 2.5, 1), fliplr(prctile(y_boot, 97.5, 1))], ...
    [0.3 0.3 0.3], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(x_fit, mean(y_boot, 1, 'omitnan'), 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
errorbar(med1, med2, mad2, mad2, mad1, mad1, 'Color', [0.4 0.4 0.4], ...
    'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.9);
scatter(med1, med2, 30, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 1);
set(gca, 'FontSize', 12);
xlabel(['log input % ' cellNames{corrPair(1)}]);
ylabel(['log input % ' cellNames{corrPair(2)}]);
axis equal;
tit = title(sprintf('r=%.2f [%.2f, %.2f], p=%.3g', r, ci_r(1), ci_r(2), p));
tit.FontWeight = 'normal';
tit.FontSize = 12;

figure;
histogram(boot_r, 50, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6]);
hold on;
xline(r, 'r-', 'LineWidth', 2, 'Label', sprintf('r=%.2f', r));
xline(ci_r(1), 'b--', 'LineWidth', 1.5);
xline(ci_r(2), 'b--', 'LineWidth', 1.5);
xlabel('Pearson r'); ylabel('Probability');
title(sprintf('Bootstrap Distribution: r=%.2f [%.2f, %.2f]', r, ci_r(1), ci_r(2)));
set(gca, 'FontSize', 12);

fprintf('Observed r: %.3f\n', r);
fprintf('Bootstrap mean r: %.3f\n', mean(boot_r, 'omitnan'));
fprintf('Bootstrap median r: %.3f\n', median(boot_r, 'omitnan'));
fprintf('Bootstrap std: %.3f\n', std(boot_r, 'omitnan'));

%% Input bias ntsr vs grp per area
R = analyze_input_bias(cnt{1}', cnt{2}', areaNames, ...
    'TotalU', D(1).nTotal', 'TotalL', D(2).nTotal', ...
    'ColorU', colors.col(5, :), 'ColorL', colors.col(6, :), ...
    'PosPalette', colors.col_bias_med, 'NegPalette', colors.col_bias_lat, ...
    'GroupNames', cellNames(1:2));
save(fullfile(P.output, 'bias_ntsr_vs_grp_pooled.mat'), 'R');
writetable(removevars(R, {'Bootstrap_Distributions', 'Permutation_Distributions'}), ...
    fullfile(P.output, 'bias_ntsr_vs_grp_pooled.csv'));

%% Hierarchical clustering of areas per cell type and Mantel test
% Areas are clustered by the Pearson correlation of their input fractions
% across brains (distance 1 - r, average linkage).
propA = cnt{1}' ./ D(1).nTotal';                 % brains x areas, ntsr
propB = cnt{2}' ./ D(2).nTotal';                 % grp
corrA = corr(propA, 'rows', 'pairwise', 'Type', 'Pearson');
corrB = corr(propB, 'rows', 'pairwise', 'Type', 'Pearson');

clustA = plot_cluster_heatmap(corrA, areaNames, clusterCutoff, 'Wide-field neurons', [100 100 541.6 519.2]);
clustB = plot_cluster_heatmap(corrB, areaNames, clusterCutoff, 'Narrow-field neurons', [750 100 541.6 519.2]);

fprintf('\n=== WIDE-FIELD CLUSTERS ===\n');
for k = 1:max(clustA), fprintf('Cluster %d: %s\n', k, strjoin(areaNames(clustA == k), ', ')); end
fprintf('\n=== NARROW-FIELD CLUSTERS ===\n');
for k = 1:max(clustB), fprintf('Cluster %d: %s\n', k, strjoin(areaNames(clustB == k), ', ')); end

% Mantel test: correlation of the two correlation matrices (upper
% triangles), null by permuting the areas of the second matrix
nA = numel(areaNames);
upperMask = triu(true(nA), 1);
triA = corrA(upperMask);
similarity = corr(triA, corrB(upperMask));
null_corrs = zeros(nPerm, 1);
rng(seed);
for k = 1:nPerm
    idx = randperm(nA);
    permB = corrB(idx, idx);
    null_corrs(k) = corr(triA, permB(upperMask));
end
p_mantel = mean(null_corrs >= similarity);
p_mantel_two = mean(abs(null_corrs) >= abs(similarity));

fprintf('\n=== MANTEL TEST ===\n');
fprintf('Structure similarity (r): %.3f\n', similarity);
fprintf('P-value (one-sided): %.4f\n', p_mantel);
fprintf('P-value (two-sided): %.4f\n', p_mantel_two);

figure('Position', [400 400 500 400]);
histogram(null_corrs, 50, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
hold on;
xline(similarity, 'r-', 'LineWidth', 2);
xlabel('Correlation between matrices');
ylabel('Count');
title(sprintf('Mantel Test: r = %.3f, p = %.4f', similarity, p_mantel));
legend({'Null distribution', 'Observed'}, 'Location', 'best');


%% Local functions
function plot_knee(knee, kneeMean, threshRange, cols)
% Cumulative distribution with the knee reference line, and distance curve
figure('Position', [368.2 465.8 624.8 208], 'Name', 'Threshold Analysis');

subplot(1, 2, 1); hold on;
yCurve = zeros(numel(knee), 1);
for k = 1:numel(knee)
    pm = sort(mean(knee(k).pct, 2), 'descend');
    cumPct = (1:numel(pm)) / numel(pm) * 100;
    plot(pm, cumPct, '-', 'LineWidth', 2, 'Color', cols(k, :));
    [pmU, iU] = unique(pm, 'stable');
    yCurve(k) = interp1(pmU, cumPct(iU), kneeMean, 'linear', 'extrap');
end
boundary = threshRange(end);
plot([0, boundary], [100, 0], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
yRef = 100 - (100 / boundary) * kneeMean;
plot([kneeMean, kneeMean], [yRef, mean(yCurve)], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
xlabel('input % threshold'); ylabel('cumulative % of areas');
xlim([-1, 20]);
set(gca, 'FontSize', 11, 'box', 'off');

subplot(1, 2, 2); hold on;
yPeak = 0;
for k = 1:numel(knee)
    plot(threshRange, knee(k).distances, '-', 'LineWidth', 2, 'Color', cols(k, :));
    plot(knee(k).thresh, knee(k).distances(knee(k).idx), 'o', 'MarkerSize', 8, ...
        'Color', cols(k, :), 'MarkerFaceColor', cols(k, :));
    yPeak = max(yPeak, interp1(threshRange, knee(k).distances, kneeMean, 'linear', 'extrap'));
end
plot([kneeMean, kneeMean], [0, yPeak], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
xlabel('input % threshold'); ylabel('dist from ref line');
xlim([-0.1, 4.7]); ylim([0.48 0.65]);
set(gca, 'FontSize', 11, 'box', 'off');
end

function plot_grouped_percent(individuals, medVals, errVals, sortedLabels, ctColors, cellNames)
% Bars (median) with error bars (MAD) and brains as dots, per area and cell
% type; background shading by major brain region
figure('Position', [462.6 354.6 735.2 216]);
nGroups = numel(individuals);
groupwidth = nGroups / (nGroups - 0.7);
lr = size(sortedLabels, 1);

for i = 1:nGroups
    x = (1:2:lr*2) - groupwidth/2 + (2*i - 1) * groupwidth / (2*nGroups);
    draw_area_bars(x, individuals{i}, medVals(:, i), errVals(:, i), ctColors(i, :), 0.7, true);
end

region = cell2mat(sortedLabels(:, 3));
shade_brain_regions(region, max(cellfun(@(v) max(v(:)), individuals)) + 3);
label_brain_regions(region, 21);

ylim([0 23]);
ax = gca;
set(ax, 'XTick', (1:2:lr*2), 'TickLabelInterpreter', 'none', ...
    'XTickLabel', sortedLabels(:, 1), 'XTickLabelRotation', 90, 'FontSize', 10, 'FontName', 'Arial', ...
    'GridColor', [.5 .5 .5], 'GridLineStyle', '--', 'YGrid', 'on', ...
    'YTick', 0:10:100, 'YTickLabel', string(0:10:100), ...
    'XGrid', 'off', 'color', 'none', 'TickDir', 'out', ...
    'GridAlpha', 0.5, 'GridLineWidth', 0.8);
ax.XAxis.TickLength = [0 0.5];
ax.YRuler.TickLabelGapOffset = +5;
curX = xlim;
xlim([0 curX(2)]);
ylabel('% of input');
box off;

legend('', 'Location', 'east');
for j = 1:numel(cellNames)
    patch([NaN NaN NaN NaN], [NaN NaN NaN NaN], ctColors(j, :), 'FaceAlpha', 1, ...
        'EdgeColor', 'none', 'DisplayName', cellNames{j});
end
hold off;
set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0.02));
set(gcf, 'Renderer', 'painters');
end

function clusters = plot_cluster_heatmap(C, areaNames, cutoffFrac, ttl, pos)
% Correlation heatmap ordered by hierarchical clustering, with dendrogram
n = numel(areaNames);
distVec = squareform(1 - C, 'tovector');
link  = linkage(distVec, 'average');
order = optimalleaforder(link, distVec);
cutoff = cutoffFrac * max(link(:, 3));
clusters = cluster(link, 'cutoff', cutoff, 'criterion', 'distance');
bounds = [0.5; find(diff(clusters(order)) ~= 0) + 0.5; n + 0.5];

heatLeft = 0.15; heatBottom = 0.15; heatSize = 0.55; dendHeight = 0.12;
figure('Position', pos);
axes('Position', [heatLeft heatBottom heatSize heatSize]);
imagesc(C(order, order));
caxis([-1 1]);
colormap(redblue(256));
colorbar('Position', [heatLeft + heatSize + 0.02, heatBottom, 0.03, heatSize]);
hold on;
for i = 1:numel(bounds) - 1
    rectangle('Position', [bounds(i), bounds(i), bounds(i+1) - bounds(i), bounds(i+1) - bounds(i)], ...
        'EdgeColor', 'k', 'LineWidth', 2);
end
set(gca, 'XTick', 1:n, 'YTick', 1:n, 'XTickLabel', areaNames(order), ...
    'YTickLabel', areaNames(order), 'XTickLabelRotation', 90);

axes('Position', [heatLeft, heatBottom + heatSize, heatSize, dendHeight]);
dendrogram(link, 0, 'Reorder', order, 'ColorThreshold', cutoff);
hold on;
yline(cutoff, 'k--', 'LineWidth', 1);
set(gca, 'XLim', [0.5 n + 0.5], 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
box off;
title(ttl, 'FontSize', 14);
end

function cmap = redblue(n)
% Blue - white - red colormap
if nargin < 1, n = 256; end
half = floor(n / 2);
r = [linspace(0, 1, half), ones(1, n - half)]';
g = [linspace(0, 1, half), linspace(1, 0, n - half)]';
b = [ones(1, half), linspace(1, 0, n - half)]';
cmap = [r g b];
end
