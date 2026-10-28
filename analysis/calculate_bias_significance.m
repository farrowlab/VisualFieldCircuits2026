%% Upper vs lower visual field bias of inputs to sSC cell types
% Loads the preprocessed Med (upper visual field) and Lat (lower visual
% field) area tables of one cell type, selects the areas to analyse,
% determines the input-percentage threshold (knee of the cumulative
% distribution), tests each area for an upper/lower visual field bias
% (analyze_input_bias) and correlates upper vs lower input percentages.
%
% Requires the output of load_data_preprocess_inputs_upper_lower.m.
% Output (P.output): bias_<celltype>.mat (Results table R) and
%                    bias_<celltype>.csv (R without the resampling distributions)

clearvars; close all;

P = setup_paths();

%% Parameters
celltype = 'ntsr';                      % 'ntsr' or 'grp'

% Area selection. An area is kept if (each condition in Med or in Lat):
%   - it passes the inclusion criteria of the preprocessing (Threshold),
%   - its mean enrichment over baseline density (TT column 15) >= minDensity,
%   - its mean % of inputs outside the sSC (TT column 9) >= minPercent,
% and it is not one of the excluded rows.
minDensity = 1.6;
minPercent = 1.11;
excludePrefix = 'SCs';                                          % all sSC layers
excludeAreas  = ["MBsen", "MBmot", "MBsta", "ICc", "ICd", "ICe", "ACAd"];

% Knee analysis of the input-percentage threshold
threshRange = 0:0.01:20;

% Upper vs lower correlation: bootstrap over brains
nBootCorr = 10000;
seedCorr  = 42;

set(0, 'defaultAxesFontSize', 10);

%% Load data
S   = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB');
BBB = S.BBB;
colors = load(fullfile(P.resources, 'color_schemes.mat'), ...
    'col', 'col_bias_med_gray', 'col_bias_lat_gray');

S  = load(fullfile(P.processed, ['OD_' celltype '.mat']), 'OD_cell');
OD = S.OD_cell;

sides = {'Med', 'Lat'};                 % Med = upper, Lat = lower visual field
for s = 1:2
    T = load(fullfile(P.processed, [celltype '_' sides{s} '.mat']), 'TT', 'Threshold', 'tableid');
    side(s).count     = cell2mat(T.TT(:, 7));         % areas x brains
    side(s).percent   = cell2mat(T.TT(:, 9));         % % of inputs outside sSC
    side(s).density   = cell2mat(T.TT(:, 15));        % enrichment over baseline density
    side(s).threshold = T.Threshold;
    side(s).nTotal    = [OD.nTotalminusSC{T.tableid}]; % inputs outside sSC per brain
end

%% Select areas
acronyms = string(BBB(:, 1));
passThreshold = side(1).threshold ~= 0 | side(2).threshold ~= 0;
passDensity   = mean(side(1).density, 2) >= minDensity | mean(side(2).density, 2) >= minDensity;
passPercent   = mean(side(1).percent, 2) >= minPercent | mean(side(2).percent, 2) >= minPercent;
excluded      = startsWith(acronyms, excludePrefix) | ismember(acronyms, excludeAreas);

keepIdx = find(passThreshold & passDensity & passPercent & ~excluded);
fprintf('Final kept areas: %d of %d\n', numel(keepIdx), numel(acronyms));

areaNames = BBB(keepIdx, 1)';
for s = 1:2
    side(s).count   = side(s).count(keepIdx, :);
    side(s).percent = side(s).percent(keepIdx, :);
    side(s).density = side(s).density(keepIdx, :);
end

sideNames = {'medial', 'lateral'};
for s = 1:2
    fprintf('Total number of included %s %s input cells: %d\n', celltype, sideNames{s}, ...
        sum(side(s).count, 'all'));
end
for s = 1:2
    fprintf('Average cell counts per %s brain (%s): %d\n', celltype, sideNames{s}, ...
        round(mean(sum(side(s).count, 1))));
end

%% Input-percentage threshold: knee of the cumulative distribution
knee = struct();
for s = 1:2
    fprintf('\nAnalyzing: %s\n', sideNames{s});
    [knee(s).thresh, knee(s).distances, knee(s).idx] = ...
        find_cumulative_knee(side(s).percent, threshRange);
end
kneeMean = mean([knee.thresh]);
fprintf('\nKnee threshold  medial: %.2f%%  lateral: %.2f%%  mean: %.2f%%\n', ...
    knee(1).thresh, knee(2).thresh, kneeMean);
for s = 1:2
    fprintf('  %s: %d areas retained\n', sideNames{s}, sum(mean(side(s).percent, 2) >= kneeMean));
end

plot_knee_analysis(side, knee, kneeMean, threshRange, sideNames);

%% Upper vs lower bias per area
R = analyze_input_bias(side(1).count', side(2).count', areaNames, ...
    'TotalU', side(1).nTotal', 'TotalL', side(2).nTotal', ...
    'ColorU', colors.col(1, :), 'ColorL', colors.col(2, :), ...
    'PosPalette', colors.col_bias_med_gray, 'NegPalette', colors.col_bias_lat_gray);

save(fullfile(P.output, ['bias_' celltype '.mat']), 'R');
writetable(removevars(R, {'Bootstrap_Distributions', 'Permutation_Distributions'}), ...
    fullfile(P.output, ['bias_' celltype '.csv']));

%% Correlation of upper vs lower input percentages
% Per area: median and median absolute deviation of log input % across brains
pctMed = side(1).percent;
pctLat = side(2).percent;
nonzeroVals = [pctMed(pctMed > 0); pctLat(pctLat > 0)];
pseudocount = min(nonzeroVals) / 10;              % below the smallest nonzero value

logMed = log(pctMed + pseudocount);
logLat = log(pctLat + pseudocount);
medianMed = median(logMed, 2);
medianLat = median(logLat, 2);
madMed = mad(logMed, 1, 2);
madLat = mad(logLat, 1, 2);

% Bootstrap: resample brains (separately for Med and Lat) for the CI of the
% Pearson r and of the regression line
nMed = size(logMed, 2);
nLat = size(logLat, 2);
finiteMed = medianMed(isfinite(medianMed));
x_fit  = linspace(min(finiteMed), max(finiteMed), 100);
boot_r = zeros(nBootCorr, 1);
y_boot = zeros(nBootCorr, numel(x_fit));

rng(seedCorr);
for b = 1:nBootCorr
    bMed = median(logMed(:, randi(nMed, nMed, 1)), 2);
    bLat = median(logLat(:, randi(nLat, nLat, 1)), 2);
    valid = isfinite(bMed) & isfinite(bLat);
    if sum(valid) > 2
        boot_r(b) = corr(bMed(valid), bLat(valid), 'type', 'Pearson');
        y_boot(b, :) = polyval(polyfit(bMed(valid), bLat(valid), 1), x_fit);
    else
        boot_r(b) = NaN;
        y_boot(b, :) = NaN;
    end
end

ci_r = prctile(boot_r, [2.5 97.5]);
y_lower = prctile(y_boot, 2.5, 1);
y_upper = prctile(y_boot, 97.5, 1);
y_fit   = median(y_boot, 1, 'omitnan');

valid = isfinite(medianMed) & isfinite(medianLat);
[r, p] = corr(medianMed(valid), medianLat(valid), 'type', 'Pearson');

figure('Position', [329.8 326.6 482.4 344.8]);
hold on;
fill([x_fit, fliplr(x_fit)], [y_lower, fliplr(y_upper)], [0.3 0.3 0.3], ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
plot(x_fit, y_fit, 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
errorbar(medianMed, medianLat, madLat, madLat, madMed, madMed, 'Color', [0.4 0.4 0.4], ...
    'LineStyle', 'none', 'CapSize', 0, 'LineWidth', 0.9);
scatter(medianMed, medianLat, 30, [0.4 0.4 0.4], 'filled', 'MarkerFaceAlpha', 1);
set(gca, 'FontSize', 12);
xlabel('log input % up NF');
ylabel('log input % low NF');
tit = title(sprintf('r=%.2f [%.2f, %.2f], p=%.3g', r, ci_r(1), ci_r(2), p));
tit.FontWeight = 'normal';
tit.FontSize = 12;

figure;
histogram(boot_r, 50, 'Normalization', 'probability', 'FaceColor', [0.6 0.6 0.6]);
hold on;
xline(r, 'r-', 'LineWidth', 2, 'Label', sprintf('r=%.2f', r));
xline(ci_r(1), 'b--', 'LineWidth', 1.5);
xline(ci_r(2), 'b--', 'LineWidth', 1.5);
xlabel('Pearson r');
ylabel('Probability');
title(sprintf('Bootstrap Distribution: r=%.2f [%.2f, %.2f]', r, ci_r(1), ci_r(2)));
set(gca, 'FontSize', 12);

fprintf('Observed r: %.3f\n', r);
fprintf('Bootstrap mean r: %.3f\n', mean(boot_r, 'omitnan'));
fprintf('Bootstrap median r: %.3f\n', median(boot_r, 'omitnan'));
fprintf('Bootstrap std: %.3f\n', std(boot_r, 'omitnan'));


%% Local functions
function plot_knee_analysis(side, knee, kneeMean, threshRange, sideNames)
% Cumulative distribution of mean input %, knee detection and areas retained
colors = [0.2 0.4 0.8; 0.9 0.4 0.2];
red    = [0.8 0.2 0.2];
figure('Position', [100 100 1600 450], 'Name', 'Threshold Analysis');

% 1: cumulative distribution
subplot(1, 4, 1); hold on;
for s = 1:2
    pm = sort(mean(side(s).percent, 2), 'descend');
    plot(pm, (1:numel(pm)) / numel(pm) * 100, '-', 'LineWidth', 1.5, ...
        'Color', colors(s, :), 'DisplayName', sideNames{s});
end
xline(kneeMean, '--', 'Color', red, 'LineWidth', 1.2, 'DisplayName', sprintf('Threshold (%.2f%%)', kneeMean));
yline(50, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', '50% of areas');
xlabel('Input Percentage (%)'); ylabel('Cumulative % of Areas');
title('Cumulative Distribution');
legend('Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3); xlim([0, 25]);

% 2: zoom with the reference line of the knee method
subplot(1, 4, 2); hold on;
for s = 1:2
    pm = sort(mean(side(s).percent, 2), 'descend');
    plot(pm, (1:numel(pm)) / numel(pm) * 100, '-', 'LineWidth', 1.5, ...
        'Color', colors(s, :), 'DisplayName', sideNames{s});
end
plot([0, threshRange(end)], [100, 0], '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
    'DisplayName', 'Reference line');
xline(kneeMean, '--', 'Color', red, 'LineWidth', 1.2, 'DisplayName', sprintf('Knee (%.2f%%)', kneeMean));
xlabel('Input Percentage (%)'); ylabel('Cumulative % of Areas');
title('Zoomed with Reference Line');
legend('Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3); xlim([0, 5]);

% 3: distance from the reference line (peak = knee)
subplot(1, 4, 3); hold on;
for s = 1:2
    plot(threshRange, knee(s).distances, '-', 'LineWidth', 1.5, ...
        'Color', colors(s, :), 'DisplayName', sideNames{s});
    plot(knee(s).thresh, knee(s).distances(knee(s).idx), 'o', 'MarkerSize', 8, ...
        'Color', colors(s, :), 'MarkerFaceColor', colors(s, :), 'HandleVisibility', 'off');
end
xline(kneeMean, '--', 'Color', red, 'LineWidth', 1.2, 'DisplayName', sprintf('Mean (%.2f%%)', kneeMean));
xlabel('Input Percentage (%)'); ylabel('Distance from Reference Line');
title('Knee Detection (Peak = Threshold)');
legend('Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3);

% 4: number of areas retained
subplot(1, 4, 4); hold on;
for s = 1:2
    pm = mean(side(s).percent, 2);
    plot(threshRange, sum(pm >= threshRange, 1), '-', 'LineWidth', 1.5, ...
        'Color', colors(s, :), 'DisplayName', sideNames{s});
end
xline(kneeMean, '--', 'Color', red, 'LineWidth', 1.2, 'DisplayName', sprintf('Threshold (%.2f%%)', kneeMean));
xlabel('Threshold (%)'); ylabel('Number of Areas');
title('Areas Retained');
legend('Location', 'northeast', 'FontSize', 8);
grid on; set(gca, 'GridAlpha', 0.3); xlim([0, 5]);

sgtitle(sprintf('Optimal Threshold = %.2f%%', kneeMean), 'FontSize', 12);
end
