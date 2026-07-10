# MOSAIX

A suite of analysis tools for mapping cortical cell types to their functional responses, combining voltage imaging and mFISH spatial transcriptomics datasets.

## Components

Each component is self-contained and has its own README with full setup and usage instructions.

| Component | What it does |
|-----------|--------------|
| [Cell Type-Response Mapping](./Cell%20Type-Response%20Mapping/) | HCR/MERFISH-based cortical cell-type classification + optogenetic electrophysiology aggregation (the MOSAIX pipeline) |
| [GEVI Timeseries Extraction](./GEVI%20Timeseries%20Extraction/) | `VoltageGUI` — MATLAB app to process voltage-imaging movies, draw ROIs, extract time series, and run spike/PSTH/GLM analyses |
| [Gene Annotation for GEVI ROIs](./Gene%20Annotation%20for%20GEVI%20ROIs/) | `annotateROI` — MATLAB app for binary (yes/no) labeling of ROIs by fluorescence brightness |
| [HCR Image Unmixing](./HCR%20Image%20Unmixing/) | `unmixV5` — MATLAB app for spectral unmixing of multichannel fluorescence images |

More analysis components are being added. How the pieces fit together — some run in sequence, others in parallel — will be documented here as a flowchart once the set is complete.

## Getting started

Open the component you need and follow its README. The components are independent, so you only need the dependencies for the ones you use.

Bulk input data for **Cell Type-Response Mapping** is too large for GitHub and is hosted separately on Zenodo — see that component's README for download instructions.
