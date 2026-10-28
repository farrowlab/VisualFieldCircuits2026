function cfg = config_local()
%CONFIG_LOCAL Machine-specific settings. Copy this file to config_local.m
%   and edit the paths below. config_local.m is ignored by git.

% Folder with the downloaded data (default: <repo>/data)
cfg.dataRoot = '/path/to/downloaded/data';

% Only needed for mapseq_export_data.m: MAPseq pipeline output with
% subfolders M1 (with v2), M2, M3
% cfg.mapseqRawRoot = '/path/to/MapNeuro/Data';
end
