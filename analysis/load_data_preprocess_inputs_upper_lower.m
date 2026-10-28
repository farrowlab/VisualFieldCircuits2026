%% Upper vs lower visual field inputs to Ntsr1+ and Grp+ sSC neurons: preprocessing
% Loads the labelled-cell coordinates of each tracing experiment, assigns
% cells to atlas areas and computes per-area counts, percentages and
% densities. Medial sSC injections (upper visual field, 'Med') and lateral
% sSC injections (lower visual field, 'Lat') are combined separately. Areas
% that pass the inclusion criteria are flagged in Threshold.
%
% Input
%   P.brainData/<experiment>.mat  table 'Coordinates', one row per labelled cell
%   P.atlas/areaGroupImport.mat   area table BBB
% Output (P.processed)
%   Coord_<celltype>.mat          per-experiment cell coordinates
%   OD_<celltype>.mat             per-experiment area tables and totals
%   <celltype>_<Med|Lat>.mat      combined area table TT, coordinates, Threshold

clearvars; close all;

P = setup_paths();

%% Parameters
celltypes = {'ntsr', 'grp'};

% Experiment sets (see experiment_list). Non-default sets are saved with
% the set name appended to the file names.
defaultSet = struct('ntsr', 'all', 'grp', 'all');
expSet     = defaultSet;

% Areas excluded from all analyses (prefix match, subdivisions included)
areasToRemove = {'RSPv', 'LGd'};

% Inclusion criteria: mean cell count >= minMeanCells and labelled in
% >= minRepeats brains (ntsr: 7 Med, 7 Lat brains; grp: 6 Med, 5 Lat brains)
minMeanCells = 2;
minRepeats   = struct('ntsr', 5, 'grp', 4);

%% Atlas area table
S   = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB');
BBB = S.BBB;
A   = area_indices(BBB);
SCid = A.SC;                                   % saved with the output
fprintf('sSC rows: %s\n', strjoin(string(BBB(SCid, 1)), ', '));

groupLabels = {'1', '2'};                      % '1' = Med, '2' = Lat

%% Process each cell type
for c = 1:numel(celltypes)
    ct = celltypes{c};
    if strcmp(expSet.(ct), defaultSet.(ct)), tag = ''; else, tag = ['_' expSet.(ct)]; end

    E = experiment_list(ct, expSet.(ct));

    OD = struct();
    OD.exp   = E(:, 1)';
    OD.type  = E(:, 2)';
    OD.group = groupLabels(1 + strcmp(OD.type, 'Lat'));

    Coord = struct();
    Coord.exp  = OD.exp;
    Coord.type = OD.type;

    % --- Load cells and count per area ---
    for i = 1:numel(OD.exp)
        Coord.ori{1, i} = load_cell_coordinates(P.brainData, OD.exp{i}, areasToRemove);

        [T, xyz, nUnassigned] = count_cells_per_area(Coord.ori{1, i}, BBB, str2double(OD.group{i}));
        OD.data{1, i}        = T;
        OD.coordinates{1, i} = xyz;
        if nUnassigned > 0
            fprintf('  %s: %d cells not assigned to any area\n', OD.exp{i}, nUnassigned);
        end
    end

    Coord_cell = Coord;
    save(fullfile(P.processed, ['Coord_' ct tag '.mat']), 'Coord_cell');

    % --- Totals per experiment ---
    for i = 1:numel(OD.exp)
        tot = input_totals(OD.data{i}, A);
        for f = fieldnames(tot)'
            OD.(f{1}){1, i} = tot.(f{1});
        end
    end

    % --- Combine brains per injection site (Med / Lat) ---
    injTypes = unique(OD.type);
    for t = 1:numel(injTypes)
        typeName = injTypes{t};
        tableid  = find(strcmp(OD.type, typeName));
        expid    = OD.exp(tableid);

        countMat   = zeros(A.nAreas, numel(tableid));
        CoordMerge = cell(A.nAreas, numel(tableid));
        for i = 1:numel(tableid)
            countMat(:, i)   = cell2mat(OD.data{tableid(i)}.GroupCount);
            CoordMerge(:, i) = OD.data{tableid(i)}.xyz_find;
        end

        D = totals_for_area_table(OD, tableid);
        [TT, TT_title] = build_area_table(BBB, countMat, D);
        Threshold = area_threshold(TT, minMeanCells, minRepeats.(ct));
        fprintf('%s %s: %d of %d areas pass the inclusion criteria\n', ...
            ct, typeName, nnz(Threshold), A.nAreas);

        save(fullfile(P.processed, [ct '_' typeName tag '.mat']), ...
            'TT', 'CoordMerge', 'TT_title', 'tableid', 'expid', 'Threshold');
    end

    OD_cell = OD;
    save(fullfile(P.processed, ['OD_' ct tag '.mat']), 'OD_cell', 'SCid');
    fprintf('%s data processed\n', ct);
end
