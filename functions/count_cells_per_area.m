function [T, xyzAssigned, nUnassigned] = count_cells_per_area(Coordinates, BBB, groupIdx)
%COUNT_CELLS_PER_AREA Assign labelled cells to the areas of the area table.
%   [T, xyzAssigned, nUnassigned] = COUNT_CELLS_PER_AREA(Coordinates, BBB, groupIdx)
%
%   Coordinates  table of labelled cells with variables avIndex (Allen
%                annotation index) and AP_location, DV_location, ML_location
%   BBB          area table (see area_indices)
%   groupIdx     numeric group label stored in T.GroupIdx
%
%   T            one row per area of BBB, with variables acronym_all,
%                name_all, GroupIdx, GroupCount (number of cells), vol_all
%                (area volume) and xyz_find (AP/DV/ML of those cells).
%                A cell is counted in every row whose annotation indices
%                contain it.
%   xyzAssigned  AP/DV/ML of cells assigned to at least one area
%   nUnassigned  number of cells not assigned to any area

nAreas      = size(BBB, 1);
acronym_all = BBB(:, 1);
name_all    = BBB(:, 3);
avIndex     = BBB(:, 4);
vol_all     = BBB(:, 8);

idx_data = [Coordinates.avIndex{:}];
xyz_data = [Coordinates.AP_location Coordinates.DV_location Coordinates.ML_location];

isAssigned = false(height(Coordinates), 1);
GroupIdx   = repmat(groupIdx, nAreas, 1);
GroupCount = cell(nAreas, 1);
xyz_find   = cell(nAreas, 1);
for j = 1:nAreas
    inArea = ismember(idx_data, avIndex{j});
    xyz_find{j}   = xyz_data(inArea, :);
    GroupCount{j} = nnz(inArea);
    isAssigned(inArea) = true;
end

T = table(acronym_all, name_all, GroupIdx, GroupCount, vol_all, xyz_find);
xyzAssigned = xyz_data(isAssigned, :);
nUnassigned = nnz(~isAssigned);
end
