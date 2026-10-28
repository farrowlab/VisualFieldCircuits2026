function [TT, TT_title] = build_area_table(BBB, countMat, D)
%BUILD_AREA_TABLE Per-area counts, percentages and densities across brains.
%   [TT, TT_title] = BUILD_AREA_TABLE(BBB, countMat, D)
%
%   BBB       area table (see area_indices)
%   countMat  nAreas-by-nExp cell counts
%   D         struct of 1-by-nExp per-experiment denominators: nTotal,
%             nTotalminusSC, nMax, nTotalminusCTX, nSC, SCDensity,
%             denTotalminusSC, denSum
%
%   TT has one row per area; columns 7-16 hold one value per brain.
%     1-5  acronym, name, parent area, parent group, volume (from BBB)
%     6    Repetions: number of brains with at least one labelled cell
%     7    cellNumber
%     8    % of total inputs                           (cellNumber / nTotal)
%     9    % of inputs outside the sSC, used in paper  (cellNumber / nTotalminusSC)
%     10   % of the area with most cells               (cellNumber / nMax)
%     11   cell density                                (cellNumber / area volume)
%     12   % of inputs outside CTX, sCTX and sSC       (cellNumber / nTotalminusCTX)
%     13   % of sSC cells                              (cellNumber / nSC)
%     14   density relative to sSC density             (density / SCDensity)
%     15   enrichment over baseline density outside sSC (density / denTotalminusSC)
%     16   density as % of summed densities outside sSC (density / denSum)

TT_title = {'arc', 'name', 'parentalArea', 'parentalIdx', 'areaVolume', 'Repetions', ...
    'cellNumber', 'cellPerTotal', 'cellPerTotalminusSC', 'cellPerMax', ...
    'cellDensity', 'cellPerTotalminuCTX', 'cellPerTotalSC', 'cellDenNorm', ...
    'cellDenTotalminusSC', 'cellDenSumTotal'};

TT = BBB(:, [1 3 5 6 8]);
for n = 1:size(BBB, 1)
    cellNumber  = countMat(n, :);
    cellDensity = cellNumber ./ TT{n, 5};

    TT{n, 6}  = nnz(cellNumber > 0);
    TT{n, 7}  = cellNumber;
    TT{n, 8}  = cellNumber ./ D.nTotal * 100;
    TT{n, 9}  = cellNumber ./ D.nTotalminusSC * 100;
    TT{n, 10} = cellNumber ./ D.nMax * 100;
    TT{n, 11} = cellDensity;
    TT{n, 12} = cellNumber ./ D.nTotalminusCTX * 100;
    TT{n, 13} = cellNumber ./ D.nSC * 100;
    TT{n, 14} = cellDensity ./ D.SCDensity;
    TT{n, 15} = cellDensity ./ D.denTotalminusSC;
    TT{n, 16} = cellDensity ./ D.denSum * 100;
end
end
