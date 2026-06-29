# MOSAIX Pipeline

Cortical cell type classification using HCR (Hybridization Chain Reaction) imaging, MERFISH spatial transcriptomics, and optogenetic electrophysiology.

## Pipeline Overview

```
Stream 1: Cell Classification → Electrophysiology
  R/run_classification.R  →  HCR_classification.csv  →  MATLAB/Find_EPSPs.m  →  MATLAB/Analyze_all_responses.m

Stream 2: scRNA-seq → GMM Thresholding → Reference Data
  R/make_meta_seq_cells.R  →  intermediate/All_genes_bin.csv + GMM_Thresholds.csv
```

## Directory Structure

```
config/
  paths.json                 - Shared directory paths
  experiments/               - Per-experiment JSON configs (genes, slice params, img_lims)
R/
  run_classification.R       - Main entry point: loads data, runs classifier, writes CSV
  classifier_functions.R     - HierarchicalGMMBayes, BayesClassifier, etc.
  config_loader.R            - Reads JSON configs and resolves paths
  hcr_data_loader.R          - Loads Ilastik CSVs and performs spatial rotation
  fit_gmm_thresholds.R       - GMM fitting (replaces Python script)
  make_meta_seq_cells.R      - Optional one-time: meta cells → GMM → binarization (outputs shipped in intermediate/)
  prepare_reference_data.R   - Optional one-time: builds processed_abc/ reference objects from raw Allen data (not distributed)
MATLAB/
  Find_EPSPs.m               - Unified electrophysiology extraction (config-driven)
  load_slice_data.m          - Per-slice file loading helper
  Analyze_all_responses.m    - Population analysis across experiments
input/
  <experiment_folder>/       - Curated per-experiment inputs (HCR gene tables, TS files, CheRiff TIFs/masks)
  processed_abc/             - MERFISH/ABC reference data
intermediate/                - Cross-language data (All_genes_bin.csv, GMM_Thresholds.csv, etc.)
output/
  <experiment_id>/           - Per-experiment outputs (HCR_classification.csv, .mat files)
  all_table.mat              - All cells with metadata and is_responder flag + all_med_norm_amps
  responders_table.mat       - Responder cells only
  non_responders_table.mat   - Non-responder cells only
_archive/                    - Original scripts preserved for reference
```

## Input File Formats

The per-experiment input data lives under `input/<experiment_folder>/`, pointed to by `mosaix_data_dir` (`input/`) in `config/paths.json` joined with each experiment's `data_subdir`. **This data is not stored in git** — it is too large for GitHub. Download it from Zenodo and extract it into the `input/` folder before running the pipeline (see [Downloading the input data](#downloading-the-input-data) below). Only the curated subset of files the pipeline actually reads is included; the bulk raw imaging data (registered TIFs, z-stacks, etc.) is not part of the archives and is not needed to run the pipeline. The sections below describe each input file's structure so it can be recreated.

Experiment folder mapping (`data_subdir`):

| Experiment | `data_subdir` | Type |
|------------|---------------|------|
| 0416 cortical | `Cortical_0416/` | HCR only (flat) |
| 1022 voltron | `Voltron_1022/` | HCR only (flat) |
| 1205 mosaix | `MOSAIX_1205/` | HCR only (flat) |
| 0730 thalamus | `MThal_0730/Registered/` | HCR + ephys |
| 0802 thalamus | `MThal_0802/Registered/` | HCR + ephys |
| 0811 contra | `Contra_0811/Registered/` | HCR + ephys |
| 0814 contra | `Contra_0814/Registered/` | HCR + ephys |
| 1008 contra | `Contra_1008/Registered/` | HCR + ephys |

### Directory layout (per ephys experiment, under `input/`)

```
<data_subdir>/                          (e.g. input/MThal_0730/Registered/)
  HCR_classification.csv               - written by run_classification.R; read by Find_EPSPs.m
  Slice N/                              - one folder per slice (matched by '*Slice*')
    Gene expression tables/            - Ilastik HCR gene tables (one CSV per gene)
      s0N_Penk..._table.csv
      s0N_Calb1..._table.csv
      ...
    MAX_channel_3_*.tif                - CheRiff channel max-projection image
    *ilastik_masks_AllVoltronPos_consec*.tif  - Voltron soma mask (one label per cell)
    TS files/
      TS*.mat                          - voltage-imaging time series
      Opto*.mat                        - optogenetic stimulus timing
      ROI*.mat                         - ROI → cell-ID mapping
      *CheRiff_composite_table*.csv    - direct-stimulation annotation
```

For the three HCR-only (flat) experiments, `data_subdir` points directly at a folder of `s##_<Gene>..._table.csv` files (no `Slice N/` subfolders, no ephys data).

### HCR imaging data — `Slice N/Gene expression tables/s##_<Gene>...table.csv`

Ilastik object-classification exports, one CSV per (slice, gene), read by `R/hcr_data_loader.R`. These are **inputs**, not pipeline outputs — the classifier consumes them and writes `HCR_classification.csv`. (For flat experiments the CSVs sit directly in `data_subdir`.) Slice index comes from the `s##` filename prefix; gene name from the rest of the filename (per `gene_file_regex`). The relevant columns are mapped through the experiment config's `file_colnames` (defaults shown):

| Config key | Default column | Meaning |
|------------|----------------|---------|
| `cell_id`  | `CellID` | Object/cell identifier (matches ROI `roiIndices`) |
| `X`/`Y`/`Z`| `Center_X`/`Center_Y`/`Center_Z` | Object centroid coordinates (pixels) |
| `Roi_mean` | `MeanROIIntensity` | Mean fluorescence intensity within the object |
| `Neighbor_mean` | `MeanNeighborhoodIntensity` | Mean intensity in the surrounding neighborhood |
| `Label` | `UserLabel` | Ilastik binary class (`Label 1` = positive for that gene, else negative) |

To recreate: segment each registered slice in Ilastik, classify objects per gene, and export the object-feature table with the columns above (or set `file_colnames` to match your export).

### Time series — `TS files/TS*.mat`

One or more `.mat` files per slice (which to load is set by `ephys.slice_file_counts.ts_idx` / `ts_concat`). Each contains a single variable:

- **`ts`** — the voltage-imaging traces as a numeric matrix of `timepoints × cells`, **or** a cell array of such matrices (the loader applies `cell2mat`). Column order corresponds to the ROI order in the matching `ROI*.mat` file.

### Optogenetic stimulus — `TS files/Opto*.mat`

Loaded per `ephys.slice_file_counts.opto_idx`. Each contains a single variable:

- **`stim`** — a binary column vector, length = number of timepoints, with `1` at every frame where the optogenetic stimulus is ON (`0` otherwise). The loader takes `find(stim == 1)` as the stimulus onset indices, offsetting by prior-file length when time series are concatenated.

### ROI map — `TS files/ROI*.mat`

Loaded per `ephys.slice_file_counts.roi_idx`. Each contains a single variable:

- **`roiIndices`** — a vector of cell IDs (one per imaged ROI), in the same order as the columns of `ts`. These IDs are matched against the HCR `cell_id`/`CellID` values to pair each trace with its cell-type classification.

### Direct-stimulation annotation — `TS files/*CheRiff_composite_table*.csv`

A table read with `readtable`; the only column used is:

- **`UserLabel`** — per-cell flag marking cells that were directly stimulated by CheRiff (indexed by kept cell ID into `cheriff_pos`).

### CheRiff laminar density — `output/<experiment_id>/s##_CheRiff_Input_to_All_Rois.csv`

**Derived** (not raw): produced by `MATLAB/Extract_CheRiff_density.m` from the CheRiff channel TIF and Voltron mask, then read back by `load_slice_data.m`. One CSV per slice, one row per ROI:

| Column | Meaning |
|--------|---------|
| `cellID` | ROI/cell label (matches the Voltron mask labels) |
| `allLayers` | Mean CheRiff signal over the full cortical column for that soma's x-strip |
| `normAllLayers` | `allLayers` normalized to the slice max (used as `cheriff_data`) |
| `L1` / `normL1` | Mean (and slice-normalized) CheRiff signal in L1 (top 80 px below cortex top) |
| `atSoma` / `normAtSoma` | Mean (and slice-normalized) CheRiff signal in a 101 px box around the soma |

To recreate: run `Extract_CheRiff_density('config/experiments/<id>.json')`, which needs `MAX_channel_3_*.tif` and `*ilastik_masks_AllVoltronPos_consec*.tif` present in each slice folder.

## Downloading the input data

The `input/` data is hosted on Zenodo (not in this repo). Download and extract it into the `input/` folder at the pipeline root so the structure matches the layout described above (`input/<experiment_folder>/`, `input/processed_abc/`, etc.).

Everything needed to run the pipeline lives in a single archive (~19 GB): the `processed_abc/` reference data plus all per-experiment folders.

| Archive | Size | Contents |
|---------|------|----------|
| **Pipeline inputs** | ~19 GB | `processed_abc/` (MERFISH/ABC reference) plus all per-experiment folders. Everything required to run classification and electrophysiology. |

After downloading, your `input/` folder should look like:

```
input/
  processed_abc/        (MERFISH/ABC reference objects)
  Cortical_0416/
  Voltron_1022/
  MOSAIX_1205/
  MThal_0730/
  MThal_0802/
  Contra_0811/
  Contra_0814/
  Contra_1008/
  HCR annotations/
```

> **Zenodo DOI / link:** _add once the record is published._

The raw Allen Brain Cell Atlas data (`raw_abc/`, ~36 GB) used by `R/prepare_reference_data.R` to build `processed_abc/` is **not distributed** — `processed_abc/` is provided directly. `prepare_reference_data.R` is included for transparency (it documents how the reference was generated); to re-run it from scratch you would need to download the raw Allen data yourself.

## Quick Start

### Stream 2: Prepare reference data (optional — outputs already shipped)
The binarized reference files this produces are committed under `intermediate/`, so you do **not** need to run this to use the pipeline. It is provided for transparency / to regenerate the reference from `processed_abc/`.
```r
# In RStudio, set working directory to the MOSAIX pipeline root
source('R/make_meta_seq_cells.R')
# Outputs (already provided in intermediate/): All_genes_bin.csv, GMM_Thresholds.csv,
#   meta_log_counts.csv, meta_norm_counts.csv
```

### Stream 1: Classify cells for an experiment
```r
# Edit experiment_config_file in R/run_classification.R, then:
source('R/run_classification.R')
# Or from command line:
# Rscript R/run_classification.R config/experiments/2024-07-30_thalamus.json
```

### Stream 1: Extract electrophysiology data
```matlab
% In MATLAB, cd to the MOSAIX pipeline root
Find_EPSPs('config/experiments/2024-07-30_thalamus.json')
% Then aggregate:
Analyze_all_responses
% Outputs: output/all_table.mat, responders_table.mat, non_responders_table.mat
```

### Working with analysis outputs

`all_table.mat` contains every cell across experiments with an `is_responder` column (1 = responder, 0 = non-responder). It also includes `all_med_norm_amps`, a vector of median-normalized response amplitudes aligned to the rows of `all_table`. Use the `is_responder` column to select responder or non-responder subsets of `all_med_norm_amps` for downstream analysis:

```matlab
load('output/all_table.mat');  % loads all_table and all_med_norm_amps
responder_amps = all_med_norm_amps(logical(all_table.is_responder));
non_responder_amps = all_med_norm_amps(~logical(all_table.is_responder));
```

## Adding a New Experiment

1. Create a JSON config in `config/experiments/` (copy an existing one as template)
2. Set the experiment-specific parameters: `data_subdir`, `slice_offset`, `n_slices`, `img_lims`, and `ephys` section
3. Run classification: `source('R/run_classification.R')` with the new config
4. Run electrophysiology: `Find_EPSPs('config/experiments/your_new_experiment.json')`

## Dependencies

**R**: dplyr, Seurat, ggplot2, jsonlite, reticulate, tidyr, hdf5r, Polychrome, igraph, reshape2, pheatmap, patchwork, cowplot, scales

**Python** (called from R via reticulate): scikit-learn, numpy

**MATLAB**: Signal Processing Toolbox (for `medfilt1`), Statistics and Machine Learning Toolbox (for `tinv`)

**Data**: Downloaded separately from Zenodo into `input/` — curated per-experiment inputs in `input/<experiment_folder>/` and reference files in `input/processed_abc/`. See [Downloading the input data](#downloading-the-input-data). No external drive needed to run the pipeline.
