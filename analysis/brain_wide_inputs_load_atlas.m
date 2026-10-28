%% Brain-wide inputs: load atlases
% Loads the Allen CCF annotation volume, the fUSi atlas and the custom area
% grouping used for the brain-wide input analysis.
%
% Requirements: see README.md (MATLAB version, toolboxes, external code).

clearvars; close all;

P = setup_paths();

%% 1. Allen CCF atlas
% av: annotation volume (int16); st: structure tree table
[av, st] = get_Allen_Data();

%% 2. fUSi atlas
% fusiAtlas.mat was generated once from subAtlas.mat with:
%   D = decomposition2(fullfile(P.atlas, 'list_areas_shortened_264_CL02.txt'), subAtlas.sa);
subAtlas = load(fullfile(P.atlas, 'subAtlas.mat'));
load(fullfile(P.atlas, 'fusiAtlas.mat'), 'D');
load(fullfile(P.atlas, 'LinMap.mat'));           % TODO: list loaded variable(s)

%% 3. Custom area grouping and plotting colours
load(fullfile(P.atlas, 'areaGroupImport.mat'));        % TODO: list loaded variable(s)
load(fullfile(P.resources, 'color_schemes.mat'));      % TODO: list loaded variable(s)
