function keepIdx = select_areas(D, minDensity, minPercent, excluded)
%SELECT_AREAS Areas included in the per-area analyses and figures.
%   keepIdx = SELECT_AREAS(D, minDensity, minPercent, excluded)
%
%   D         struct array, one element per cell type and/or injection site,
%             with fields threshold (inclusion flag from the preprocessing),
%             density (areas x brains, TT column 15) and percent (TT column 9)
%   excluded  logical mask of rows that are never kept
%
%   An area is kept if, in at least one element of D, it passes the
%   inclusion criteria (threshold ~= 0); in at least one, its mean density
%   is >= minDensity; and in at least one, its mean percent is >= minPercent.

n = numel(D(1).threshold);
passThr = false(n, 1);
passDen = false(n, 1);
passPct = false(n, 1);
for k = 1:numel(D)
    passThr = passThr | D(k).threshold ~= 0;
    passDen = passDen | mean(D(k).density, 2) >= minDensity;
    passPct = passPct | mean(D(k).percent, 2) >= minPercent;
end
keepIdx = find(passThr & passDen & passPct & ~excluded(:));
end
