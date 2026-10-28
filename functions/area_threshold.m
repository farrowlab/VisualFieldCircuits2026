function Threshold = area_threshold(TT, minMeanCells, minRepeats)
%AREA_THRESHOLD Flag areas that pass the inclusion criteria.
%   Threshold = AREA_THRESHOLD(TT, minMeanCells, minRepeats) returns 1 for
%   areas whose mean cell count across brains (TT column 7) is at least
%   minMeanCells and that are labelled in at least minRepeats brains
%   (TT column 6), and 0 otherwise.

meanCells = cellfun(@mean, TT(:, 7));
repeats   = [TT{:, 6}]';
Threshold = double(meanCells >= minMeanCells & repeats >= minRepeats);
end
