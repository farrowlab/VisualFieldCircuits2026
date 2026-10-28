%% Ipsilateral vs contralateral inputs to Ntsr1+ and Grp+ sSC neurons: preprocessing
% Same as load_data_preprocess_inputs_upper_lower.m, but the labelled cells
% of each brain are split into the hemisphere ipsilateral (ML >= 0) and
% contralateral (ML < 0) to the injection. Percentages are relative to all
% inputs of the brain (both hemispheres) minus the ipsilateral sSC, so the
% contralateral sSC counts as an input area.
%
% Input
%   P.brainData/<experiment>.mat  table 'Coordinates', one row per labelled cell
%   P.atlas/areaGroupImport.mat   area table BBB
% Output (P.processed)
%   Coord_<celltype>_ipsicontra.mat            cell coordinates per brain and hemisphere
%   OD_<celltype>_ipsicontra.mat               per-brain/hemisphere area tables and totals
%   <celltype>_<Med|Lat>_<ipsi|contra>.mat     combined area table TT, coordinates, Threshold
%
% Entries of OD and Coord.ori alternate per brain:
%   brain 1 ipsi, brain 1 contra, brain 2 ipsi, brain 2 contra, ...

clearvars; close all;

P = setup_paths();

%% Parameters
celltypes = {'ntsr', 'grp'};

% Experiment sets (see experiment_list). For ntsr only brains whose
% injection-site labelling is confined to the ipsilateral sSC are used.
% Non-default sets are saved with the set name appended to the file names.
defaultSet = struct('ntsr', 'ipsiOnly', 'grp', 'all');
expSet     = defaultSet;

% Areas excluded from all analyses (prefix match, subdivisions included)
areasToRemove = {'RSPv', 'LGd'};

% Inclusion criteria: mean cell count >= minMeanCells and labelled in
% >= minRepeats brains (ntsr: 5 Med, 7 Lat brains; grp: 6 Med, 5 Lat brains)
minMeanCells = 2;
minRepeats   = struct('ntsr', 4, 'grp', 4);

%% Atlas area table
S   = load(fullfile(P.atlas, 'areaGroupImport.mat'), 'BBB');
BBB = S.BBB;
A   = area_indices(BBB);
SCid = A.SC;                                   % saved with the output
fprintf('sSC rows: %s\n', strjoin(string(BBB(SCid, 1)), ', '));

% Rows used for the summed density (excludes sSC and MB summary rows)
notMBorSC = A.notMB;
notMBorSC(SCid) = false;

% Annotation indices of the sSC, for the hemisphere sanity check
scAnnotation = cell2mat(cellfun(@(x) x(:), BBB(SCid, 4), 'UniformOutput', false));

sides       = {'ipsi', 'contra'};
groupLabels = {'1', '2'};                      % '1' = Med, '2' = Lat
groupSuffix = {'', '.5'};                      % ipsi: '1', '2'; contra: '1.5', '2.5'

%% Process each cell type
for c = 1:numel(celltypes)
    ct = celltypes{c};
    if strcmp(expSet.(ct), defaultSet.(ct)), tag = ''; else, tag = ['_' expSet.(ct)]; end

    E       = experiment_list(ct, expSet.(ct));
    nBrains = size(E, 1);
    nEntries = 2 * nBrains;

    OD    = struct();
    Coord = struct();
    Coord.exp  = E(:, 1)';                     % one name per brain
    Coord.type = E(:, 2)';

    % --- Load cells, split by hemisphere, count per area ---
    for b = 1:nBrains
        C = load_cell_coordinates(P.brainData, E{b, 1}, areasToRemove);
        isIpsi = C.ML_location >= 0;

        % Sanity check: the injected (ipsilateral) hemisphere should hold
        % most sSC cells; a low value means the brain is mirrored in ML
        sSCcells = ismember([C.avIndex{:}], scAnnotation);
        fprintf('  %s: %.0f%% of sSC cells ipsilateral\n', E{b, 1}, ...
            100 * nnz(sSCcells(:) & isIpsi) / nnz(sSCcells));

        baseGroup = groupLabels{1 + strcmp(E{b, 2}, 'Lat')};
        for s = 1:2
            k = 2 * (b - 1) + s;
            if s == 1, Cs = C(isIpsi, :); else, Cs = C(~isIpsi, :); end

            Coord.ori{1, k} = Cs;
            OD.exp{1, k}    = [E{b, 1} '_' sides{s}];
            OD.type{1, k}   = [E{b, 2} '_' sides{s}];
            OD.group{1, k}  = [baseGroup groupSuffix{s}];

            [T, xyz, nUnassigned] = count_cells_per_area(Cs, BBB, str2double(OD.group{k}));
            OD.data{1, k}        = T;
            OD.coordinates{1, k} = xyz;
            if nUnassigned > 0
                fprintf('  %s: %d cells not assigned to any area\n', OD.exp{k}, nUnassigned);
            end
        end
    end

    save(fullfile(P.processed, ['Coord_' ct '_ipsicontra' tag '.mat']), 'Coord');

    % --- Totals per brain and hemisphere ---
    nTotalNoIC = zeros(1, nEntries);           % used for the whole-brain total
    for k = 1:nEntries
        count = cell2mat(OD.data{k}.GroupCount);
        vol   = cell2mat(OD.data{k}.vol_all);

        nTotalNoIC(k)          = sum(count(A.notMBorIC));
        OD.nTotal{1, k}        = sum(count(A.notMB));
        OD.nSC{1, k}           = sum(count(SCid));
        OD.nCTX{1, k}          = sum(count(A.CTX));
        OD.nsCTX{1, k}         = sum(count(A.sCTX));
        OD.nTotalminusCTX{1, k} = OD.nTotal{k} - OD.nCTX{k} - OD.nsCTX{k} - OD.nSC{k};
        OD.nTH{1, k}           = sum(count(A.TH));
        OD.nMax{1, k}          = max(count);
        OD.aSC{1, k}           = sum(vol(SCid));
        OD.tot_vol{1, k}       = sum(vol(A.notMB));
        OD.aTotalminusSC{1, k} = OD.tot_vol{k} - OD.aSC{k};
    end

    % --- Whole-brain totals (both hemispheres), stored for both entries ---
    for b = 1:nBrains
        kIpsi   = 2 * b - 1;
        kContra = 2 * b;
        nTOTAL     = nTotalNoIC(kIpsi) + nTotalNoIC(kContra);
        nSC_ipsi   = OD.nSC{kIpsi};
        nSC_contra = OD.nSC{kContra};
        for k = [kIpsi kContra]
            OD.nTOTAL{1, k}        = nTOTAL;
            OD.nSC_ipsi{1, k}      = nSC_ipsi;
            OD.nSC_contra{1, k}    = nSC_contra;
            OD.nSCtot{1, k}        = nSC_ipsi + nSC_contra;
            OD.nTotalminusSC{1, k} = nTOTAL - nSC_ipsi;     % all inputs minus ipsilateral sSC
            OD.SCDensity{1, k}     = nSC_ipsi / OD.aSC{1};  % sSC volume is the same in all brains
        end
    end

    % --- Densities ---
    for k = 1:nEntries
        count = cell2mat(OD.data{k}.GroupCount);
        vol   = cell2mat(OD.data{k}.vol_all);
        OD.denTotalminusSC{1, k} = OD.nTotalminusSC{k} / OD.aTotalminusSC{k};
        OD.denSum{1, k} = sum(count(notMBorSC) ./ vol(notMBorSC), 'omitnan');
    end

    % --- Combine brains per injection site and hemisphere ---
    entryTypes = unique(OD.type);
    for t = 1:numel(entryTypes)
        typeName = entryTypes{t};
        tableid  = find(strcmp(OD.type, typeName));
        expid    = OD.exp(tableid);

        countMat   = zeros(A.nAreas, numel(tableid));
        CoordMerge = cell(A.nAreas, numel(tableid));
        for i = 1:numel(tableid)
            countMat(:, i)   = cell2mat(OD.data{tableid(i)}.GroupCount);
            CoordMerge(:, i) = OD.data{tableid(i)}.xyz_find;
        end

        D.nTotal          = [OD.nTOTAL{tableid}];
        D.nTotalminusSC   = [OD.nTotalminusSC{tableid}];
        D.nMax            = [OD.nMax{tableid}];
        D.nTotalminusCTX  = [OD.nTotalminusCTX{tableid}];
        D.nSC             = [OD.nSC_ipsi{tableid}];
        D.SCDensity       = [OD.SCDensity{tableid}];
        D.denTotalminusSC = [OD.denTotalminusSC{tableid}];
        D.denSum          = [OD.denSum{tableid}];

        [TT, TT_title] = build_area_table(BBB, countMat, D);
        Threshold = area_threshold(TT, minMeanCells, minRepeats.(ct));
        fprintf('%s %s: %d of %d areas pass the inclusion criteria\n', ...
            ct, typeName, nnz(Threshold), A.nAreas);

        save(fullfile(P.processed, [ct '_' typeName tag '.mat']), ...
            'TT', 'CoordMerge', 'TT_title', 'tableid', 'expid', 'Threshold');
    end

    save(fullfile(P.processed, ['OD_' ct '_ipsicontra' tag '.mat']), 'OD', 'SCid');
    fprintf('%s data processed\n', ct);
end
