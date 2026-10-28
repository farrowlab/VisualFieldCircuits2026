%% MAPseq: export barcode matrices from the sequencing pipeline output
% Builds the data files shared with the paper from the MAPseq pipeline
% output of the three adult mice (M1-M3). Needs the pipeline output, which
% is not public; it documents how the shared data were made.
%
% Per mouse, the pipeline output folder contains a .mat file (barcode
% matrix B<m>, spike-normalised matrix B<m>norm, spike-in reads 'spikes')
% and an .xlsx sample sheet (column of each sample in the matrix, sample
% names). Mouse 1 was processed in subfolder v2 and its samples are
% renamed and reordered to match mice 2 and 3.
%
% Output: P.mapseq/mapseq_barcode_matrices.xlsx
%   sheets M1norm, M2norm, M3norm   spike-normalised counts, one row per
%                                   barcode, one column per sample (all
%                                   barcodes, before any filtering)
%   sheet samples                   per sample: name in each mouse, number
%                                   of spike-in molecules, and the factor
%                                   that converts normalised values back to
%                                   molecule counts (count = norm * scale)

clearvars; close all;

P = setup_paths();
assert(isfield(P, 'mapseqRaw') && isfolder(P.mapseqRaw), ...
    'Set cfg.mapseqRawRoot in config_local.m to the MAPseq pipeline output folder.');

nSamples = 39;            % samples analysed (mouse 1 has extra samples, dropped)
outFile  = fullfile(P.mapseq, 'mapseq_barcode_matrices.xlsx');
if ~isfolder(P.mapseq), mkdir(P.mapseq); end

%% Load pipeline output of each mouse
B = cell(1, 3); Bnorm = cell(1, 3); spikes = cell(1, 3);
sampleNum = cell(1, 3); sampleName = cell(1, 3);
for m = 1:3
    folder = fullfile(P.mapseqRaw, sprintf('M%d', m));
    if m == 1, folder = fullfile(folder, 'v2'); end

    matFile = dir(fullfile(folder, '*.mat'));
    fprintf('M%d: %s\n', m, matFile(1).name);
    S = load(fullfile(folder, matFile(1).name));
    B{m}     = S.(sprintf('B%d', m));
    Bnorm{m} = S.(sprintf('B%dnorm', m));
    spikes{m} = cellfun(@length, {S.spikes.counts2u});   % spike-in molecules per sample
    others = setdiff(intersect(fieldnames(S), {'B1', 'B2', 'B3'}), sprintf('B%d', m));
    if ~isempty(others)
        warning('%s also contains %s; only B%d is used.', matFile(1).name, strjoin(others, ', '), m);
    end

    % Sample sheet: column of each sample in the matrix, and sample names
    sheet = dir(fullfile(folder, '*xlsx'));
    T = readtable(fullfile(folder, sheet(1).name));
    hdr = T.Properties.VariableNames;
    colVar = hdr(ismember(hdr, {'Column_InMatrix', 'columnsInMatrix'}));
    sampleNum{m} = T.(colVar{1});
    if ismember('Brain', hdr)
        sampleNum{m} = sampleNum{m}(T.Brain == m);
    end
    names = table2array(T(:, 2));
    if m > 1, names = names(2:end); end      % first row of the sheets of mice 2 and 3 is not a sample
    sampleName{m} = names(1:numel(sampleNum{m}));
end

%% Mouse 1: rename samples and reorder them to match mice 2 and 3
sampleName{1}([8 9 18 19 24 37]) = {'SCsci', 'SCscc', 'SCdci', 'SCdcc', 'SCdami', 'SCdaalc'};
order = (1:numel(sampleName{1}))';
for i = 1:nSamples
    order(i) = find(strcmp(sampleName{1}, sampleName{2}{i}));
end
sampleName{1} = sampleName{1}(order);
sampleNum{1}  = sampleNum{1}(order);

%% Sample columns in the order of the sample sheets
for m = 1:3
    cols = sampleNum{m}(1:nSamples);
    B{m}      = B{m}(:, cols);
    Bnorm{m}  = Bnorm{m}(:, cols);
    spikes{m} = spikes{m}(cols);
    sampleName{m} = sampleName{m}(1:nSamples);
end

if ~isequal(sampleName{1}, sampleName{2}, sampleName{3})
    warning('Sample names differ between mice; check that the columns match:');
    disp(table((1:nSamples)', sampleName{1}, sampleName{2}, sampleName{3}, ...
        'VariableNames', {'column', 'M1', 'M2', 'M3'}));
end

%% Factor converting normalised values back to molecule counts
scale = cell(1, 3);
for m = 1:3
    scale{m} = norm_to_count_scale(B{m}, Bnorm{m});
    err = max(abs(round(Bnorm{m} .* scale{m}) - B{m}), [], 'all');
    fprintf('M%d: %d barcodes; counts recovered from normalised values with max error %g\n', ...
        m, size(B{m}, 1), err);
    if err > 0
        warning('M%d: normalised values do not convert exactly back to counts.', m);
    end
end

%% Write
if isfile(outFile), delete(outFile); end
for m = 1:3
    writetable(array2table(Bnorm{m}, 'VariableNames', sampleName{m}'), outFile, ...
        'Sheet', sprintf('M%dnorm', m));
end
samples = table((1:nSamples)', sampleName{1}, sampleName{2}, sampleName{3}, ...
    spikes{1}', spikes{2}', spikes{3}', scale{1}', scale{2}', scale{3}', ...
    'VariableNames', {'column', 'name_M1', 'name_M2', 'name_M3', ...
    'spikein_M1', 'spikein_M2', 'spikein_M3', 'scale_M1', 'scale_M2', 'scale_M3'});
writetable(samples, outFile, 'Sheet', 'samples');
fprintf('Written: %s\n', outFile);


%% Local functions
function scale = norm_to_count_scale(B, Bnorm)
% Per-sample factor with B = Bnorm .* scale (median ratio over nonzero entries)
scale = zeros(1, size(B, 2));
for j = 1:size(B, 2)
    nz = Bnorm(:, j) ~= 0;
    if any(nz)
        scale(j) = median(B(nz, j) ./ Bnorm(nz, j));
    end
end
end
