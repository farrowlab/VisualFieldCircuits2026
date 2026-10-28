%% Upper vs lower visual field bias of inputs: Ntsr1+ vs Grp+ sSC neurons
% Loads the preprocessed Med (upper visual field) and Lat (lower visual
% field) area tables of both cell types, selects one common set of areas,
% and tests per area the upper/lower bias of each cell type and the
% difference between cell types (analyze_input_bias_two_celltypes).
%
% Requires the output of load_data_preprocess_inputs_upper_lower.m.
% Output (P.output): bias_ntsr_vs_grp.mat (Results table R) and .csv

clearvars; close all;

P = setup_paths();

%% Parameters
celltypes = {'ntsr', 'grp'};            % CT1, CT2

% Area selection, combined across cell types and injection sites. An area
% is kept if, in at least one cell type and site, it passes the inclusion
% criteria of the preprocessing (Threshold); in at least one, its mean
% enrichment (TT column 15) >= minDensity; in at least one, its mean % of
% inputs outside the sSC (TT column 9) >= minPercent; and it is not one of
% the excluded rows.
minDensity = 1.6;
minPercent = 1.11;
excludePrefix = 'SCs';                                          % all sSC layers
excludeAreas  = ["MBsen", "MBmot", "MBsta", "ICc", "ICd", "ICe", "ACAd"];

set(0, 'defaultAxesFontSize', 10);

%% Load data
S   = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB');
BBB = S.BBB;
colors = load(fullfile(P.resources, 'color_schemes.mat'), 'col');

sides = {'Med', 'Lat'};                 % Med = upper, Lat = lower visual field
for c = 1:numel(celltypes)
    S  = load(fullfile(P.processed, ['OD_' celltypes{c} '.mat']), 'OD_cell');
    OD = S.OD_cell;
    for s = 1:2
        T = load(fullfile(P.processed, [celltypes{c} '_' sides{s} '.mat']), 'TT', 'Threshold', 'tableid');
        D(c, s).count     = cell2mat(T.TT(:, 7));         % areas x brains
        D(c, s).percent   = cell2mat(T.TT(:, 9));         % % of inputs outside sSC
        D(c, s).density   = cell2mat(T.TT(:, 15));        % enrichment over baseline density
        D(c, s).threshold = T.Threshold;
        D(c, s).nTotal    = [OD.nTotalminusSC{T.tableid}]; % inputs outside sSC per brain
    end
end

%% Select areas (same set for both cell types)
acronyms = string(BBB(:, 1));
passThreshold = false(numel(acronyms), 1);
passDensity   = false(numel(acronyms), 1);
passPercent   = false(numel(acronyms), 1);
for c = 1:numel(celltypes)
    for s = 1:2
        passThreshold = passThreshold | D(c, s).threshold ~= 0;
        passDensity   = passDensity   | mean(D(c, s).density, 2) >= minDensity;
        passPercent   = passPercent   | mean(D(c, s).percent, 2) >= minPercent;
    end
end
excluded = startsWith(acronyms, excludePrefix) | ismember(acronyms, excludeAreas);

keepIdx = find(passThreshold & passDensity & passPercent & ~excluded);
fprintf('Final kept areas: %d of %d\n', numel(keepIdx), numel(acronyms));
areaNames = BBB(keepIdx, 1)';

for c = 1:numel(celltypes)
    for s = 1:2
        D(c, s).count = D(c, s).count(keepIdx, :);
    end
    fprintf('Total number of included %s medial input cells: %d\n',  celltypes{c}, sum(D(c, 1).count, 'all'));
    fprintf('Total number of included %s lateral input cells: %d\n', celltypes{c}, sum(D(c, 2).count, 'all'));
    fprintf('Average cell counts per %s brain (medial): %d\n',  celltypes{c}, round(mean(sum(D(c, 1).count, 1))));
    fprintf('Average cell counts per %s brain (lateral): %d\n', celltypes{c}, round(mean(sum(D(c, 2).count, 1))));
end

%% Upper vs lower bias per cell type and cell-type difference
R = analyze_input_bias_two_celltypes( ...
    D(1, 1).count', D(1, 2).count', D(2, 1).count', D(2, 2).count', areaNames, ...
    'TotalU1', D(1, 1).nTotal', 'TotalL1', D(1, 2).nTotal', ...
    'TotalU2', D(2, 1).nTotal', 'TotalL2', D(2, 2).nTotal', ...
    'ColorU1', colors.col(3, :), 'ColorL1', colors.col(4, :), ...
    'ColorU2', colors.col(1, :), 'ColorL2', colors.col(2, :), ...
    'ColorIntPos', colors.col(5, :), 'ColorIntNeg', colors.col(6, :));

save(fullfile(P.output, 'bias_ntsr_vs_grp.mat'), 'R');
writetable(R, fullfile(P.output, 'bias_ntsr_vs_grp.csv'));
