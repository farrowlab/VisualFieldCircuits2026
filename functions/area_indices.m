function A = area_indices(BBB)
%AREA_INDICES Rows of the area table used to compute input totals.
%   A = AREA_INDICES(BBB) returns row indices and masks for the area table
%   BBB (loaded from atlas/areaGroupImport.mat).
%
%   Columns of BBB: 1 acronym, 3 name, 4 Allen annotation indices,
%   5 parent area, 6 parent group (1 cortex, 2 cortical subplate,
%   4 thalamus), 8 volume.
%
%   A.SC          superficial SC (injection site)
%   A.CTX, A.sCTX, A.TH   cortex, cortical subplate, thalamus
%   A.MBsen, A.MBmot, A.MBsta, A.IC   midbrain summary rows and inferior
%                 colliculus subdivisions, excluded from input totals
%   A.notMB       mask: all rows except MBsen, MBmot, MBsta
%   A.notMBorIC   mask: all rows except MBsen, MBmot, MBsta, ICc, ICd, ICe

A.nAreas = size(BBB, 1);
acronyms = string(BBB(:, 1));
assert(numel(unique(acronyms)) == A.nAreas, 'area_indices:duplicateAcronyms', ...
    'Area acronyms in BBB are not unique.');

A.SC = find(contains(string(BBB(:, 3)), 'Superior colliculus superficial layers'))';

parentGroup = [BBB{:, 6}];
A.CTX  = find(parentGroup == 1);
A.sCTX = find(parentGroup == 2);
A.TH   = find(parentGroup == 4);

A.MBsen = find(contains(acronyms, "MBsen"));
A.MBmot = find(contains(acronyms, "MBmot"));
A.MBsta = find(contains(acronyms, "MBsta"));
A.IC    = find(contains(acronyms, ["ICc", "ICd", "ICe"]));

A.notMB = true(A.nAreas, 1);
A.notMB([A.MBsen; A.MBmot; A.MBsta]) = false;
A.notMBorIC = A.notMB;
A.notMBorIC(A.IC) = false;
end
