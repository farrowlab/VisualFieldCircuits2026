%% Inputs to Ntsr1+, Grp+ and Gad2+ sSC neurons, all injections pooled: preprocessing
% Same as load_data_preprocess_inputs_upper_lower.m, but without separating
% upper and lower visual field: all brains of a cell type are combined.
% Cell coordinates per brain are read from one combined file per dataset
% instead of from the individual experiment files.
%
% Input
%   P.brainData/Coord_ntsr_grp_gad2_cvs.mat  variable Coord: 1-by-3 cell
%                                            (ntsr, grp, gad2), each with
%                                            fields exp, type, ori (cell tables)
%   P.atlas/areaGroupImport.mat              area table BBB
% Output (P.processed)
%   Coord_<celltype>_pooled.mat   per-brain cell coordinates
%   OD_<celltype>_pooled.mat      per-brain area tables and totals
%   <celltype>_pooled.mat         combined area table TT, coordinates, Threshold

clearvars; close all;

P = setup_paths();

%% Parameters
celltypes = {'ntsr', 'grp', 'gad2'};           % order of the cells in Coord

coordFile = 'Coord_ntsr_grp_gad2_cvs.mat';

% Areas excluded from all analyses (prefix match, subdivisions included)
areasToRemove = {'RSPv', 'LGd', 'TEa'};

% Inclusion criteria: mean cell count >= minMeanCells and labelled in
% >= minRepeats brains (ntsr: 8 brains, grp: 7 brains, gad2: 4 brains)
minMeanCells = 2;
minRepeats   = struct('ntsr', 5, 'grp', 5, 'gad2', 4);

%% Atlas area table
S   = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB');
BBB = S.BBB;
A   = area_indices(BBB);
SCid = A.SC;                                   % saved with the output
fprintf('sSC rows: %s\n', strjoin(string(BBB(SCid, 1)), ', '));

%% Cell coordinates of all brains
S = load(fullfile(P.brainData, coordFile), 'Coord');
CoordAll = S.Coord;
assert(numel(CoordAll) == numel(celltypes), ...
    '%s holds %d cell types, expected %d.', coordFile, numel(CoordAll), numel(celltypes));

%% Process each cell type
for c = 1:numel(celltypes)
    ct = celltypes{c};

    Coord = CoordAll{c};
    OD = struct();
    OD.exp  = Coord.exp;
    OD.type = Coord.type;
    fprintf('%s: %s\n', ct, strjoin(string(OD.exp(:))', ', '));   % check the cell-type order

    % --- Remove excluded areas and count cells per area ---
    for i = 1:numel(Coord.ori)
        Coord.ori{i} = remove_areas(Coord.ori{i}, areasToRemove);

        [T, xyz, nUnassigned] = count_cells_per_area(Coord.ori{i}, BBB, 1);
        OD.data{1, i}        = T;
        OD.coordinates{1, i} = xyz;
        if nUnassigned > 0
            fprintf('  %s: %d cells not assigned to any area\n', OD.exp{i}, nUnassigned);
        end
    end

    Coord_cell = Coord;
    save(fullfile(P.processed, ['Coord_' ct '_pooled.mat']), 'Coord_cell');

    % --- Totals per brain ---
    for i = 1:numel(OD.data)
        tot = input_totals(OD.data{i}, A);
        for f = fieldnames(tot)'
            OD.(f{1}){1, i} = tot.(f{1});
        end
    end

    % --- Combine all brains ---
    tableid = 1:numel(OD.data);
    expid   = OD.exp(tableid);

    countMat   = zeros(A.nAreas, numel(tableid));
    CoordMerge = cell(A.nAreas, numel(tableid));
    for i = 1:numel(tableid)
        countMat(:, i)   = cell2mat(OD.data{tableid(i)}.GroupCount);
        CoordMerge(:, i) = OD.data{tableid(i)}.xyz_find;
    end

    D = totals_for_area_table(OD, tableid);
    [TT, TT_title] = build_area_table(BBB, countMat, D);
    Threshold = area_threshold(TT, minMeanCells, minRepeats.(ct));
    fprintf('%s: %d of %d areas pass the inclusion criteria\n', ct, nnz(Threshold), A.nAreas);

    save(fullfile(P.processed, [ct '_pooled.mat']), ...
        'TT', 'CoordMerge', 'TT_title', 'tableid', 'expid', 'Threshold');

    OD_cell = OD;
    save(fullfile(P.processed, ['OD_' ct '_pooled.mat']), 'OD_cell', 'SCid');
    fprintf('%s data processed\n', ct);
end
