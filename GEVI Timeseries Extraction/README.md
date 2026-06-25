# VoltageGUI

`VoltageGUI.mlapp` is a MATLAB **App Designer** application for processing and
analyzing voltage-imaging movies. It is the interactive front end of the
SpikeTriggeredPlots pipeline: you load raw imaging movies, draw/refine regions
of interest (ROIs) over individual neurons, extract their fluorescence time
series, and then visualize spikes, stimulus-aligned responses, and population
analyses (PSTH, GLM, cross-correlation, etc.).

It wraps the command-line `Pipeline_*` functions in this repo behind a two-tab
GUI so the whole motion-correction → ROI → time-series → plotting workflow can
be driven by clicking, with results cached to disk so a session can be re-opened
without recomputing.

---

## Requirements

- **MATLAB R2018a or newer** (the app declares `minimumSupportedMATLABRelease =
  R2018a`). A recent release (R2021b+) is recommended since the code uses
  `width`/`height` on images, `clim`, and name=value calls.
- Toolboxes used by the underlying pipeline: **Image Processing Toolbox**,
  **Signal Processing Toolbox**, and **Statistics and Machine Learning Toolbox**
  (the GLM/PSTH/bootstrap code).
- The **entire repository on the MATLAB path.** The app calls
  `addpath(genpath(pwd))` at startup, so it must be launched **from the repo
  root** (the folder containing `VoltageGUI.mlapp`, `support/`, `Pipeline_*.m`,
  etc.). It also tries to add a sibling `../Microscope Control` folder for
  `selectROI2`; that line fails silently if the folder is absent.

### Supported input data

The file picker accepts `.tif` / `.tiff`, `.raw`, and `.prd` movies (or a folder
of numbered `*.tiff` frames). Frames are read through `FrameReader.m`. Optional
WaveSurfer ephys/stim files are auto-associated when present.

---

## Launching

From the repository root, either:

```matlab
cd '/path/to/SpikeTriggeredPlots-main'
VoltageGUI            % opens the app
```

or double-click `VoltageGUI.mlapp` in the MATLAB **Current Folder** browser
(make sure the Current Folder is the repo root first).

On startup the app immediately prompts you to **select an input file** (the
`Select New File` flow runs automatically).

---

## What it produces / where things are saved

For each selected movie, results are cached in an **`Analysis/`** subfolder next
to the movie (created automatically). File names encode the source movie, e.g.
`ROI(<movie>)_1-1.mat`, `TS(<movie>)_1-1.mat`. Key outputs:

| File token | Contents |
|------------|----------|
| `ExpInfo(...)`  | Experiment metadata generated on file load |
| `Movavg(...).tif` | Median/average reference image for the main display |
| `corrim`, `dfRangeim`, `dfRatioim`, `bgDotProdImage` | Alternate background images for ROI drawing |
| `RegBGMov(...).tif` | Motion-corrected, background-subtracted movie |
| `ROI(...).mat`  | ROI masks, per-ROI time series, neuropil range, region |
| `TS(...).mat`   | Extracted fluorescence time series |
| `OptoData(...)` | Detected light-artifact / stimulus timing |
| `GLM(...)`, `BootstrapInfo(...)` | Analysis outputs |

Because outputs are cached, the app skips recomputation when nothing relevant
has changed (tracked internally via a `hasChanged` flag). Changing motion
correction, background-ball radius, neuropil range, or the ROIs marks the data
as changed and forces re-extraction on the next plot.

---

## The interface

The window has two tabs.

### Tab 1 — Image / ROI editing

This is where you pick data, choose a reference image, and define ROIs.

**File / dataset controls**
- **Select New File** — open the file picker and add one or more movies. All
  loaded movies appear in the **File** dropdown; switching the dropdown loads
  that dataset.
- **Remove Current** — drop the currently selected dataset from the list.
- **Background Image dropdown** — choose what the main axes shows as the
  backdrop for drawing ROIs: mean image, correlation image, df range, df ratio,
  or background dot-product image. (Computed on demand and cached.)
- **Color Scheme** — colormap used for ROI contours/labels.
- **White / Black** — contrast (percentile) limits for the displayed image.
- **lock daspect** — lock the image aspect ratio to 1:1.

**Motion correction / preprocessing**
- **Motion (checkbox)** — enable rigid registration for this dataset.
- **BG ball radius** — rolling-ball background-subtraction radius (0 = off).
- **Motion Threshold**, **Temporal Filter**, **neuropil dist (min/max)** —
  parameters for extraction and motion regression.
- **Run Motion on All** — runs motion correction + background removal **and**
  time-series extraction across *every* loaded dataset, with a progress bar.
  (Requires the Motion checkbox to be on.) Will reuse the current ROIs for other
  movies of the same frame size.

**ROI creation / editing** (operate on the main image)
- **Manual ROI** + the **draw-shape dropdown** — draw one or more ROIs by hand
  (polygon/freehand/etc.). **Manual Add** merges a hand-drawn region into the
  selected ROI.
- **Smart ROI** (toggle) — click a pixel and the app grows an ROI automatically
  around it; keep clicking to add more, toggle off when done. **Smart Add**
  merges a smart region into the selected ROI.
- **Optimize spx** — optimizes ROI pixel selection to improve the signal.
- **Try Unify ROIs** — copies the current dataset's ROIs onto all other loaded
  datasets (same frame size) and re-extracts their time series. **Overwrites**
  other datasets' ROIs — it asks for confirmation first.
- **ROI List** — selectable list of ROIs; selecting one highlights its contour
  and shows its trace in the small time-series axis.
- **Delete ROIs** — remove the selected ROIs.
- **Save ROIs** / **Load ROIs** — write/read the ROI `.mat` file.
- **Frames disregard** — number of leading frames to ignore.

**Plot / PSTH**
- **Plot** — extract time series (if needed) and draw the full traces on Tab 2.
- **PSTH** — extract time series, detect spikes (using **Spike STD** threshold),
  and produce peri-stimulus / spike-triggered average plots in separate figures.

### Tab 2 — Time-series / analysis

Shows the extracted traces and analysis plots.

- **Big Time Series / Motion axes** — the main trace display (linked x-axis).
- **smooth / Spike STD / Seperation** — display smoothing, spike-detection
  threshold (in SDs), and vertical spacing between stacked traces.
- **flip dff** — invert the ΔF/F sign (voltage indicators often go negative on
  depolarization).
- **HP freq** — high-pass cutoff applied before plotting/spike detection.
- **Motion Regression dropdown** — choose whether/how motion is regressed out of
  the traces.
- **Detect light artifact** — re-detects stimulus/opto light artifacts and
  rewrites the `OptoData` file.
- **Split Stim** (toggle) — split/average traces by stimulus condition; reveals
  the stim-shift buttons for reordering stimulus correspondence.
- **lock plot** — reuse the cached plot instead of recomputing when nothing has
  changed.
- **popout** — move the trace axis into its own resizable figure.
- **GLM Plot** — fit and plot a generalized linear model of spiking for the
  selected neurons.
- **Open ROI stats**, **Open Stim Plots** — launch the companion apps
  `ROIStats.mlapp` and `time_series_plots.mlapp` on the loaded datasets.

---

## Typical workflow

1. **Launch** `VoltageGUI` from the repo root; pick one or more movies in the
   prompt.
2. Choose a **Background Image** that best reveals neurons (often correlation or
   df-ratio).
3. (Optional) Turn on **Motion** and set **BG ball radius**, then **Run Motion
   on All** if you have several recordings.
4. Draw neurons with **Smart ROI** and/or **Manual ROI**; tidy up with
   **Delete**, **Manual/Smart Add**, and **Optimize spx**.
5. Click **Save ROIs**.
6. Click **Plot** (or **PSTH**) — the app extracts time series via
   `Pipeline_getTimeSeries` and shows the traces on Tab 2.
7. On Tab 2, tune **Spike STD**, **HP freq**, **smooth**, motion regression, and
   stimulus splitting; then run **GLM Plot** / **PSTH** for analyses.
8. To process a whole set identically, define ROIs on one recording and use
   **Try Unify ROIs** or **Run Motion on All** to propagate and batch-extract.

---

## Related files in this repo

- `Pipeline_getTimeSeries.m`, `Pipeline_MotionCorrectRemoveBackground.m`,
  `Pipeline_getExpInfo.m` — the pipeline steps the GUI invokes.
- `simple_plot.m`, `spike_avg_plots.m`, `stim_avg_plot.m`, `GLMPlot` —
  plotting/analysis backends behind **Plot**, **PSTH**, and **GLM Plot**.
- `FrameReader.m`, `Select_input_files.m`, `support/getFileName.m` — data I/O
  and the `Analysis/` file-naming convention.
- `ROIStats.mlapp`, `time_series_plots.mlapp` — companion apps opened from the
  GUI.

---

## Notes & known limitations

- The app **must** be started from the repo root, since it relies on
  `genpath(pwd)` to find its dependencies.
- A known issue noted in the source: when loading ROIs saved in another session,
  the small preview time-series for the image tab may not recalculate and can
  show a stale trace until you re-extract (Plot/PSTH).
- Heavy steps (motion correction, time-series extraction over many files) show a
  wait bar and can take a while; results are cached so subsequent runs are fast
  unless inputs change.
