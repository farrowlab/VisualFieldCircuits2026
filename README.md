# VisualFieldCircuits2026

MATLAB analysis code for brain-wide input mapping and MAPseq analyses.

## Setup

1. Clone this repository.
2. In MATLAB, run `setup_paths` from the repository root. This adds `functions/`, `atlas/`, `resources/` and `external/` to the path and creates `output/` and `data/processedTable/` if needed.
3. (Optional) Copy `config_local_example.m` to `config_local.m` and set `cfg.dataRoot` if the data live outside `<repo>/data`. `config_local.m` is git-ignored.

## Data

Per-experiment `.mat` files, coordinate tables and `allen_structure_id.json` are included in `data/`.

## Atlas files (not in the repository)

The following files exceed GitHub's file-size limit and must be downloaded separately into `atlas/`:

- `annotation_volume_10um_by_index.npy`
- `annotation_volume_10um_merged_cortical_areas.npy`
- `template_volume_10um.npy`
- `annotation_10.nrrd`, `annotation_100.nrrd`

TODO: add download link (e.g. Allen CCF 2017 volumes / Zenodo).

## Layout

- `analysis/` – analysis and plotting scripts (incl. `analysis/mapseq/`)
- `functions/` – helper functions and bundled toolboxes (gramm, effect_size_matlab, matplotlib colormaps)
- `atlas/` – Allen CCF annotation tables and small atlas resources
- `external/` – third-party code (allenCCF, MousePFC_ResearchSynthesis, Registration_Segmentation_KF)
- `resources/` – shared resources such as `color_schemes.mat`
