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
  make_meta_seq_cells.R      - Full Stream 2: meta cells → GMM → binarization
  prepare_reference_data.R   - One-time: creates Seurat reference objects from Allen data
MATLAB/
  Find_EPSPs.m               - Unified electrophysiology extraction (config-driven)
  load_slice_data.m          - Per-slice file loading helper
  Analyze_all_responses.m    - Population analysis across experiments
intermediate/                - Cross-language data (All_genes_bin.csv, GMM_Thresholds.csv, etc.)
output/
  <experiment_id>/           - Per-experiment outputs (HCR_classification.csv, .mat files)
  all_table.mat              - All cells with metadata and is_responder flag + all_med_norm_amps
  responders_table.mat       - Responder cells only
  non_responders_table.mat   - Non-responder cells only
_archive/                    - Original scripts preserved for reference
```

## Quick Start

### Stream 2: Prepare reference data (run once)
```r
# In RStudio, set working directory to the MOSAIX pipeline root
source('R/make_meta_seq_cells.R')
# Outputs: intermediate/All_genes_bin.csv, intermediate/GMM_Thresholds.csv
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

**Data**: Requires access to `/Volumes/SynMap Data/` (external drive) and reference files in `input/processed_abc/`
