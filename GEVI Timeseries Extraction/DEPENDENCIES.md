# VoltageGUI — Self-contained bundle

This folder contains `VoltageGUI.mlapp`, its `README.md`, and every in-repo file
the app needs to run. It was produced by walking the app's call graph (including
the code embedded inside the `.mlapp` files) and copying each reachable
dependency while preserving the original directory layout.

## How to run

Launch **from this folder** (the app calls `addpath(genpath(pwd))` at startup,
so the whole tree must be the current folder / on the path):

```matlab
cd '.../VoltageGUI_standalone'
VoltageGUI
```

See `README.md` for full usage instructions.

## What's included

- **`VoltageGUI.mlapp`** — the app itself.
- **Companion apps** opened from the GUI: `ROIStats.mlapp`,
  `time_series_plots.mlapp`, `selectROI2.mlapp`, `getROICorrespondence.mlapp`.
- **Pipeline / processing functions** used by the app:
  `Pipeline_getExpInfo`, `Pipeline_getTimeSeries`,
  `Pipeline_MotionCorrectRemoveBackground`, `FrameReader`, `Select_input_files`,
  `rollingBall`, `getCorrelationim`, `getSTDRatioRange`, `detectLightArtifact`,
  the stim/ROI loaders (`load_stim_series`, `load_stim_correspondence`,
  `load_stim_distances`, `load_time_series`), the plotting backends
  (`simple_plot`, `spike_avg_plots`), TIFF I/O (`Fast_BigTiff_Write`,
  `saveastiff`), `brewermap`, `uipickfiles`, and others.
- **`support/`** — the helper functions the above depend on (`getFileName`,
  `getAnalysisDir`, `getAllTimeSeries`, smoothing, masks, etc.), plus:
  - **`support/GLMFiles/`** — the GLM analysis backend (`GLMPlot` and friends).
  - **`support/+ScanImageTiffReader/`** — MATLAB package used by `FrameReader`
    to read ScanImage TIFFs (copied whole, including its `private/` folder).
  - **`support/Wavesurfer-master/+ws/`** — the WaveSurfer MATLAB package, needed
    for the ephys/stimulus path (`load_stim_series` calls `ws.loadDataFile`).
    Only the `+ws` package is included, not the rest of the WaveSurfer
    distribution. This is the bulk of the folder size (~14 MB); it is only
    exercised when you load WaveSurfer stim/ephys data.

## Not included (intentionally)

Parts of the original repo unrelated to running VoltageGUI were left out — e.g.
the wide-field (`Pipeline_*WF*`) pipeline, alternate `Pipeline_master*` scripts,
standalone cross-correlation/Granger plotting scripts, `annotateROI.mlapp`, and
editor autosave (`.asv`) / macOS metadata files.

## Known caveat

The app calls a helper named **`saveTiff`** (used when caching the reference
image and df images), but no `saveTiff.m` exists anywhere in the original
repository — it appears to be an external function expected on the user's MATLAB
path. It is therefore **not** in this bundle. The closely related `saveastiff.m`
is included; if you hit an "Undefined function `saveTiff`" error, provide a
`saveTiff` function (e.g. a thin wrapper over `saveastiff`) on the path.
