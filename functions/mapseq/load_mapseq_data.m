function D = load_mapseq_data(xlsxFile)
%LOAD_MAPSEQ_DATA Read the shared MAPseq barcode matrices of mice M1-M3.
%   D = LOAD_MAPSEQ_DATA(xlsxFile) reads the file written by
%   mapseq_export_data.m and returns, per mouse m (cells of length 3):
%     D.Bnorm{m}   spike-normalised counts, barcodes x samples
%     D.B{m}       molecule counts (Bnorm .* scale, rounded)
%     D.spikes{m}  spike-in molecules per sample (1 x samples)
%     D.names{m}   sample names (samples x 1 cell)
%
%   Molecule counts are needed because the soma threshold and the
%   detection threshold in the targets are defined in molecules. They are
%   recovered from the normalised values with the per-sample factors in
%   the 'samples' sheet.

sheets = sheetnames(xlsxFile);
D = struct('Bnorm', {cell(1, 3)}, 'B', {cell(1, 3)}, 'spikes', {cell(1, 3)}, 'names', {cell(1, 3)});
for m = 1:3
    T = readtable(xlsxFile, 'Sheet', sprintf('M%dnorm', m), 'VariableNamingRule', 'preserve');
    D.Bnorm{m} = table2array(T);
    D.names{m} = T.Properties.VariableNames(:);
end

if ~ismember('samples', sheets)
    error('load_mapseq_data:noScale', ['Sheet ''samples'' not found in %s. It holds the ' ...
        'per-sample factors that convert normalised values to molecule counts, which the ' ...
        'soma and target thresholds are defined in.'], xlsxFile);
end
S = readtable(xlsxFile, 'Sheet', 'samples', 'VariableNamingRule', 'preserve');
for m = 1:3
    D.spikes{m} = S.(sprintf('spikein_M%d', m))';
    D.B{m} = round(D.Bnorm{m} .* S.(sprintf('scale_M%d', m))');
end

if ~isequal(D.names{1}, D.names{2}, D.names{3})
    warning('load_mapseq_data:names', 'Sample names differ between mice.');
end
end
