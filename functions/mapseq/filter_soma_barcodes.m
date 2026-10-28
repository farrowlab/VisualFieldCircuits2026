function [B, Bnorm, rows] = filter_soma_barcodes(B, Bnorm, somaCols, somaThresh, ratioThresh, excludeCols, originalIndexing)
%FILTER_SOMA_BARCODES Keep barcodes with a clear soma in the given samples.
%   [B, Bnorm, rows] = FILTER_SOMA_BARCODES(B, Bnorm, somaCols, somaThresh,
%                      ratioThresh, excludeCols, originalIndexing)
%
%   B, Bnorm     barcodes x samples, molecule counts and normalised values
%   somaCols     samples a soma may be in (here the ipsilateral sSC samples)
%   somaThresh   minimum molecules in the soma sample
%   ratioThresh  maximum ratio of the strongest other sample to the soma
%   excludeCols  samples ignored when looking for the strongest other sample;
%                [] = all samples except the soma
%
%   The soma of a barcode is the sample with most molecules. A barcode is
%   kept if its soma is in somaCols with at least somaThresh molecules and
%   (strongest other sample) / soma < ratioThresh. rows are the kept rows.
%
%   originalIndexing = true reproduces the first version of the analysis
%   code, in which positions in the list of barcodes passing somaThresh were
%   used as row numbers of the full matrix. This selects the intended
%   barcodes only if those passing somaThresh are the first rows of B.

[somaVal, somaCol] = max(B, [], 2);
if originalIndexing
    m = somaCol(somaVal >= somaThresh);
    cand = find(ismember(m, somaCols));
else
    cand = find(somaVal >= somaThresh & ismember(somaCol, somaCols));
end

data = B(cand, :);
[mx, mxIdx] = max(data, [], 2);
others = data;
if isempty(excludeCols)
    others(sub2ind(size(others), (1:size(others, 1))', mxIdx)) = 0;
else
    others(:, excludeCols) = 0;
end
ratio = max(others, [], 2) ./ mx;

rows = cand(ratio < ratioThresh);
B = B(rows, :);
Bnorm = Bnorm(rows, :);
end
