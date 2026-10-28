%% Labelled input cells on atlas slices
% Plots the labelled input cells of one cell type on Allen CCF slices.
% Slices can be coronal, horizontal or sagittal; cells within +-halfWidth
% of each slice are shown.
%   colorBy = 'side':  medial (upper visual field) and lateral (lower visual
%                      field) injections on the same slices, one colour each
%   colorBy = 'brain': one figure per injection site, one colour per brain
%
% Requires the output of load_data_preprocess_inputs_upper_lower.m, the
% Allen CCF files (get_Allen_Data), allenCCF (allenCCFbregma) and
% prettyslice.m on the path.

clearvars; close all;

P = setup_paths();

%% Parameters
celltype  = 'ntsr';
plane     = 'coronal';            % 'coronal', 'horizontal' or 'sagittal'
% Slice positions (mm): coronal = AP from bregma (negative = posterior);
% horizontal = DV depth from the brain surface; sagittal = ML from midline
positions = [0 -1 -2 -3 -4 -5];
halfWidth = 0.05;                 % mm: cells within +-50 um of the slice
colorBy   = 'side';               % 'side' or 'brain'

areas          = [];              % restrict to these areas, e.g. ["VISp", "VISl"]; [] = all
excludeAreas   = ["SCs", "MBsen", "MBmot", "MBsta", ...    % sSC row, midbrain summary rows
                  "MH", "LH", "CL"];                        % likely inputs to intermediate/deep SC starters
applyThreshold = true;            % leave out areas failing the inclusion criteria
oneHemisphere  = false;           % keep only cells with ML >= 0

jitterSD = 1 / 500;               % mm: random jitter so overlapping cells stay visible
seed     = 42;

%% Atlas and cell coordinates
[av, st] = get_Allen_Data();      % annotation volume (10 um voxels) and structure tree
bregma = allenCCFbregma();        % bregma in voxels: [AP DV ML]

colors = load(fullfile(P.resources, 'color_schemes.mat'), 'col');
sideColors  = colors.col(1:2, :);
% ColorBrewer Set1 (without yellow and grey), one colour per brain
brainColors = [228 26 28; 55 126 184; 77 175 74; 152 78 163; 255 127 0; 166 86 40; 247 129 191] / 255;

sides = {'Med', 'Lat'};           % Med = upper, Lat = lower visual field
cells = cell(1, 2);               % cells{side}{brain} = n x [AP DV ML] (mm)
for s = 1:2
    T = load(fullfile(P.processed, [celltype '_' sides{s} '.mat']), 'TT', 'Threshold', 'CoordMerge');
    acr  = string(T.TT(:, 1));
    keep = ~ismember(acr, excludeAreas);
    if ~isempty(areas),  keep = keep & ismember(acr, areas);   end
    if applyThreshold,   keep = keep & T.Threshold ~= 0;       end
    for i = 1:size(T.CoordMerge, 2)
        cells{s}{i} = cell2mat(T.CoordMerge(keep, i));
    end
end

%% Plot
rng(seed);
nSlices = numel(positions);
layout  = [2, ceil(nSlices / 2)];

switch colorBy
    case 'side'
        figure;
        tiledlayout(layout(1), layout(2), 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:nSlices
            ax = nexttile;
            draw_slice(av, st, bregma, plane, positions(k), ax);
            for s = 1:2
                for i = 1:numel(cells{s})
                    [x, y] = cells_in_slice(cells{s}{i}, plane, positions(k), halfWidth, ...
                        bregma, jitterSD, oneHemisphere);
                    sc = scatter(x, y, 5, 'filled', 'CData', sideColors(s, :), 'MarkerFaceAlpha', 0.7);
                    if s == 2, uistack(sc, 'down', 2); end
                end
            end
            finish_axes(ax, plane, positions(k));
        end
        set(gcf, 'renderer', 'Painters', 'WindowState', 'maximized');

    case 'brain'
        for s = 1:2
            figure;
            tiledlayout(layout(1), layout(2), 'TileSpacing', 'compact', 'Padding', 'compact');
            for k = 1:nSlices
                ax = nexttile;
                draw_slice(av, st, bregma, plane, positions(k), ax);
                for i = 1:numel(cells{s})
                    [x, y] = cells_in_slice(cells{s}{i}, plane, positions(k), halfWidth, ...
                        bregma, jitterSD, oneHemisphere);
                    c = brainColors(mod(i - 1, size(brainColors, 1)) + 1, :);
                    sc = scatter(x, y, 10, 'filled', 'CData', c, 'MarkerFaceAlpha', 0.6);
                    if s == 2, uistack(sc, 'down', 2); end
                end
                finish_axes(ax, plane, positions(k));
            end
            set(gcf, 'renderer', 'Painters', 'WindowState', 'maximized');
        end

    otherwise
        error('colorBy must be ''side'' or ''brain''.');
end


%% Local functions
function draw_slice(av, st, bregma, plane, pos, ax)
% Atlas slice at pos (mm); image rows/columns as used by cells_in_slice
switch plane
    case 'coronal'      % DV x ML
        img = squeeze(av(round(bregma(1) - pos * 100), :, :));
    case 'horizontal'   % AP x ML
        img = squeeze(av(:, round(bregma(2) + pos * 100), :));
    case 'sagittal'     % DV x AP
        img = squeeze(av(:, :, round(bregma(3) - pos * 100)))';
    otherwise
        error('plane must be ''coronal'', ''horizontal'' or ''sagittal''.');
end
prettyslice(img, st, ax);
hold(ax, 'on');
end

function [x, y] = cells_in_slice(xyz, plane, pos, halfWidth, bregma, jitterSD, oneHemisphere)
% Cells within +-halfWidth of the slice, in image (voxel) coordinates.
% xyz: n x [AP DV ML] in mm (AP positive = posterior, as in the preprocessing)
xyz = xyz + randn(size(xyz)) * jitterSD;
AP = xyz(:, 1); DV = xyz(:, 2); ML = xyz(:, 3);
switch plane
    case 'coronal',    sel = abs(AP - (-pos)) <= halfWidth;
    case 'horizontal', sel = abs(DV - pos) <= halfWidth;
    case 'sagittal',   sel = abs(ML - pos) <= halfWidth;
end
if oneHemisphere, sel = sel & ML >= 0; end

switch plane
    case 'coronal'
        x = bregma(3) - ML(sel) * 100;  y = bregma(2) + DV(sel) * 100;
    case 'horizontal'
        x = bregma(3) - ML(sel) * 100;  y = bregma(1) + AP(sel) * 100;
    case 'sagittal'
        x = bregma(1) + AP(sel) * 100;  y = bregma(2) + DV(sel) * 100;
end
end

function finish_axes(ax, plane, pos)
names = struct('coronal', 'AP', 'horizontal', 'DV', 'sagittal', 'ML');
set(ax, 'Box', 'off', 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
xlabel(ax, sprintf('%s %.2f', names.(plane), pos), 'Color', 'k');
end
