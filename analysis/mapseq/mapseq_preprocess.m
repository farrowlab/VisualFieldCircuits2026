%% MAPseq: select sSC neurons and group target samples
% Reads the barcode matrices of mice M1-M3, keeps barcodes (neurons) with
% a clear soma in the ipsilateral superficial SC (sSC), and adds grouped
% dSC and contralateral sSC target columns.
%
% Input:  P.mapseq/mapseq_barcode_matrices.xlsx (see mapseq_export_data.m)
% Output: P.processed/mapseq_preprocessed.mat, used by mapseq_analysis.m
%
% Matrices (rows = barcodes, columns = samples 1-39 + grouped columns 40-54)
%   B1, B2, B3            molecule counts per mouse
%   B1norm, B2norm, B3norm  spike-normalised per mouse
%   BigB                  counts of all mice (rows: M1, then M2, then M3)
%   BigBnorm              BigB divided by the spike-in molecules summed over
%                         the three mice (one factor per sample)
% Soma sample of each barcode: m1, m2, m3 (per mouse), maxes (BigB rows)

clearvars; close all;

P = setup_paths();

%% Parameters
% Neuron selection, step 1: soma in the ipsilateral sSC with >= somaThresh
% molecules, and the strongest other sample < ratioThreshAll x soma
somaThresh     = 30;
ratioThreshAll = 0.85;
% Step 2: strongest sample outside the ipsilateral sSC and dSC < 1/10 of soma
ratioThreshTargets = 0.1;
sSCipsiCols = [2 3 5 6 7 9 10 11 14];
dSCipsiCols = [16 17 20 21 24 25 28];

% true: reproduce the row indexing of the first version of this code
% (see filter_soma_barcodes); false: corrected indexing
originalIndexing = false;

% Grouped target columns appended as columns 40-54
groups = { ...
    'dSCam',  [17 21]; 'dSCm',  25; 'dSCal',  [16 20]; 'dSCl',  24; 'dSCc',  28; ...   % ipsilateral dSC
    'dSCamc', [18 22]; 'dSCmc', 26; 'dSCalc', [19 23]; 'dSClc', 27; 'dSCcc', 29; ...   % contralateral dSC
    'sSCamc', [4 8];   'sSCmc', 12; 'sSCalc', [5 9];   'sSClc', 13; 'sSCcc', 15};      % contralateral sSC

% Long-distance targets: sample column -> label
targetLabels = {30, 'MLR'; 32, 'Pbg'; 34, 'LP'; 36, 'dLGN'; 38, 'vLGN'};

%% Load
D = load_mapseq_data(fullfile(P.mapseq, 'mapseq_barcode_matrices.xlsx'));
SampleName = D.names;

% Ipsilateral sSC samples: name contains 'SCs' and ends with 'i' (column vector)
sSC = find(contains(SampleName{1}, 'SCs') & endsWith(SampleName{1}, 'i'));
if ~isequal(sSC(:)', sSCipsiCols)
    warning('sSC samples found by name (%s) differ from sSCipsiCols (%s).', ...
        mat2str(sSC(:)'), mat2str(sSCipsiCols));
end

% All mice together
spikesAdults = D.spikes{1} + D.spikes{2} + D.spikes{3};
BigB     = [D.B{1}; D.B{2}; D.B{3}];
BigBnorm = BigB ./ spikesAdults;
mouseIDall = [ones(size(D.B{1}, 1), 1); 2 * ones(size(D.B{2}, 1), 1); 3 * ones(size(D.B{3}, 1), 1)];

%% Select sSC neurons
% Per mouse, and for all mice together (same criteria, applied row-wise)
Bm = D.B; Bmn = D.Bnorm;
for m = 1:3
    [Bm{m}, Bmn{m}, r1] = filter_soma_barcodes(Bm{m}, Bmn{m}, sSC, somaThresh, ratioThreshAll, [], originalIndexing);
    [Bm{m}, Bmn{m}, r2] = filter_soma_barcodes(Bm{m}, Bmn{m}, sSC, somaThresh, ratioThreshTargets, ...
        unique([sSCipsiCols dSCipsiCols]), originalIndexing);
    fprintf('M%d: %d of %d barcodes kept\n', m, numel(r2), size(D.B{m}, 1));
end
[BigB, BigBnorm, r1] = filter_soma_barcodes(BigB, BigBnorm, sSC, somaThresh, ratioThreshAll, [], originalIndexing);
[BigB, BigBnorm, r2] = filter_soma_barcodes(BigB, BigBnorm, sSC, somaThresh, ratioThreshTargets, ...
    unique([sSCipsiCols dSCipsiCols]), originalIndexing);
mouseID = mouseIDall(r1(r2));
fprintf('All mice: %d barcodes kept\n', size(BigB, 1));

if ~originalIndexing
    assert(isequal(BigB, [Bm{1}; Bm{2}; Bm{3}]), 'Selection differs between per-mouse and pooled matrices.');
end

% Check against the original indexing
for m = 1:3
    [~, ~, a1] = filter_soma_barcodes(D.B{m}, D.Bnorm{m}, sSC, somaThresh, ratioThreshAll, [], true);
    Bo = D.B{m}(a1, :);
    [~, ~, a2] = filter_soma_barcodes(Bo, Bo, sSC, somaThresh, ratioThreshTargets, ...
        unique([sSCipsiCols dSCipsiCols]), true);
    [~, ~, c1] = filter_soma_barcodes(D.B{m}, D.Bnorm{m}, sSC, somaThresh, ratioThreshAll, [], false);
    Bc = D.B{m}(c1, :);
    [~, ~, c2] = filter_soma_barcodes(Bc, Bc, sSC, somaThresh, ratioThreshTargets, ...
        unique([sSCipsiCols dSCipsiCols]), false);
    if isequal(a1(a2), c1(c2))
        fprintf('M%d: original and corrected indexing select the same %d barcodes\n', m, numel(c2));
    else
        fprintf('M%d: original and corrected indexing select DIFFERENT barcodes (%d vs %d, %d in common)\n', ...
            m, numel(a2), numel(c2), numel(intersect(a1(a2), c1(c2))));
    end
end

B1 = Bm{1}; B2 = Bm{2}; B3 = Bm{3};
B1norm = Bmn{1}; B2norm = Bmn{2}; B3norm = Bmn{3};

% Soma sample of each kept barcode (row vectors)
[~, m1] = max(B1, [], 2);       m1 = m1';
[~, m2] = max(B2, [], 2);       m2 = m2';
[~, m3] = max(B3, [], 2);       m3 = m3';
[~, maxes] = max(BigB, [], 2);  maxes = maxes';

%% Grouped target columns (40-54)
add = @(X) [X, cell2mat(cellfun(@(c) sum(X(:, c), 2), groups(:, 2)', 'UniformOutput', false))];
B1 = add(B1); B2 = add(B2); B3 = add(B3);
B1norm = add(B1norm); B2norm = add(B2norm); B3norm = add(B3norm);
BigB = add(BigB); BigBnorm = add(BigBnorm);

%% Labels and target sets
sample_names = [SampleName{3}', groups(:, 1)'];
sample_names([targetLabels{:, 1}]) = targetLabels(:, 2)';

targetSets = struct( ...
    'longIpsi',       [34 32 36 38 30], ...                  % LP, Pbg, dLGN, vLGN, MLR
    'longContra',     [35 33 37 39 31], ...
    'longIpsiContra', [34 32 36 38 30 35 33 37 39 31], ...
    'dSCipsi',        40:44, ...
    'dSCcontra',      45:49, ...
    'dSCipsiContra',  40:49, ...
    'longIpsi_dSCipsi', [34 32 36 38 30 40:44]);
orderOfSC = 1:7;

%% Save
settings = struct('somaThresh', somaThresh, 'ratioThreshAll', ratioThreshAll, ...
    'ratioThreshTargets', ratioThreshTargets, 'originalIndexing', originalIndexing);
save(fullfile(P.processed, 'mapseq_preprocessed.mat'), ...
    'B1', 'B2', 'B3', 'B1norm', 'B2norm', 'B3norm', 'BigB', 'BigBnorm', 'mouseID', ...
    'm1', 'm2', 'm3', 'maxes', 'sSC', 'SampleName', 'sample_names', 'targetSets', ...
    'orderOfSC', 'settings');
fprintf('Saved %s\n', fullfile(P.processed, 'mapseq_preprocessed.mat'));
