# annotateROI

`annotateROI.mlapp` is a MATLAB **App Designer** application for **assigning a
binary label (yes / no) to each ROI** in an ROI volume, based on how bright that
ROI is in an accompanying fluorescent image relative to its local surroundings.

It is built for voltage/calcium-imaging ROI curation: you load a labeled ROI
volume together with a fluorescent reference image, and for every ROI the app
computes an intensity metric (ROI signal vs. surrounding neighborhood), sorts the
ROIs by that metric, and presents a zoomed montage of each ROI so you can quickly
mark it as a real labeled cell (**y**) or not (**n**). The result is exported as
a CSV table of per-ROI labels.

---

## Requirements

- **MATLAB R2023a or newer recommended.** The packaged app declares a minimum of
  R2018a, but it uses several newer components/functions — most notably a
  **`RangeSlider`** (R2023a+) for contrast, plus `tiledlayout`/`nexttile`
  (R2019b+) and `clim` (R2022a+). On older releases it will not open.
- **Image Processing Toolbox** — `regionprops3`, `imdilate`, `strel`,
  `tiffreadVolume`.
- **Statistics and Machine Learning Toolbox** — `prctile`.
- The bundled **`uipickfiles.m`** (included here) must be on the MATLAB path —
  it provides the file-selection dialogs.

### Inputs

The app reads **four image volumes** (TIFF stacks or `.mat`). Multi-channel
TIFFs are detected automatically and you'll be prompted to pick a channel.

| Prompt | Property | Purpose |
|--------|----------|---------|
| **Image file 1** | `imageVol` | The fluorescent reference image used to compute each ROI's intensity metric (ROI vs. neighborhood). |
| **Voltron ROI file** | `ROIVolume` | The labeled volume of the ROIs **you are annotating** — each ROI is a distinct integer ID; 0 = background. |
| **ALL ROI file** | `ROIVolumeAll` | A labeled volume of **all** detected ROIs/segments. Used to exclude *other* ROIs from the neighborhood ring when measuring background, so a neighbor's signal doesn't contaminate the metric. |
| **voltron vol file** | `voltronVolume` | A second image volume shown **side-by-side** with image 1 in each montage tile (visual cross-reference). |

If a single dataset uses one segmentation, the "Voltron ROI file" and "ALL ROI
file" can be the same file.

---

## Launching

From the folder containing these files:

```matlab
cd '.../annotateROI_standalone'
annotateROI                 % prompts you to pick the 4 input files
```

You can also pass the four file paths to skip the dialogs:

```matlab
annotateROI(imageFile1, voltronROIFile, allROIFile, voltronVolFile)
```

---

## How the metric works

For each ROI the app takes the ROI's pixels in `imageVol` and a **neighborhood
ring** around it (the ROI mask dilated by **Neighborhood sz**, with pixels
belonging to any other ROI removed). It then compares the two with the metric
chosen in the **Metric** dropdown:

1. **Z-score** — `(prctile(ROI, ROI%) − prctile(neighborhood, Neigh%)) / (std(neighborhood)/√nROIpix)`
2. **Z-score mean** — same but using means instead of percentiles
3. **Ratio** — `prctile(ROI, ROI%) / prctile(neighborhood, Neigh%)`
4. **Ratio mean** — `mean(ROI) / mean(neighborhood)`

The **Percentile ROI** and **Neighborhood** fields set the percentiles used by
the percentile-based metrics. ROIs are listed sorted by metric value, so the most
"cell-like" ROIs cluster at one end.

Each ROI becomes a row in an internal table with columns:
`CellID, Center_X, Center_Y, Center_Z, MeanROIIntensity,
MeanNeighborhoodIntensity, MetricValue, UserLabel, roiSeen, isDoublet`.

---

## The interface

- **ROIs list box** — every ROI as `ROI # <id> (<metric>)`, sorted by metric.
  Selecting one shows it in the montage. Labeled ROIs are colored
  **green** (yes) / **red** (no); doublets get a warning icon.
- **Montage panel** — a 3×3 tiled layout: the large tile is the max-projection of
  the whole field of view; the remaining tiles show the ~5 central z-slices
  around the current ROI. In each tile the **reference image (left)** and the
  **voltron volume (right)** are shown side-by-side with the ROI outlined in red.
- **Metric** dropdown, **Percentile ROI**, **Neighborhood**, **Neighborhood sz**
  — metric type and its parameters (changing them recomputes and re-sorts).
- **FOV size** — size (px) of the crop shown around each ROI.
- **CLim slider** — display contrast (color limits) for the tiles.
- **Forward / Back** — move through the ROI list (also **→/↓** and **←/↑** keys).

### Labeling controls

| Action | Button | Key | Effect |
|--------|--------|-----|--------|
| Label as **yes** | `y` | **y** | `UserLabel = 1` (green), marks ROI seen, advances |
| Label as **no**  | `n` | **n** | `UserLabel = 0` (red), marks ROI seen, advances |
| Toggle **doublet** | toggle Doublet | **d** | flips the `isDoublet` flag (warning icon) |
| **Skip non-doublets** | checkbox | — | navigation jumps only between ROIs flagged as doublets |

### Saving / restoring

- **Export** — writes the full table to a **CSV** (`uiputfile`). Before saving it
  also **extrapolates** your manual labels to ROIs you never reviewed: any unseen
  ROI with a metric **above** the highest "yes" metric is set to `1`, and any
  with a metric **below** the lowest "no" metric is set to `0`. Reviewed ROIs
  keep your manual labels.
- **Load** — import a previously exported CSV to restore labels. You'll be asked
  to map CSV columns to the required `CellID` (and optional `UserLabel`,
  `isDoublet`) fields.

---

## Typical workflow

1. Launch `annotateROI` and select the four input volumes (reference image, the
   ROIs to label, the all-ROI segmentation, and the voltron volume).
2. Pick a **Metric** and tune **Percentile ROI / Neighborhood / Neighborhood sz**
   so the list ordering separates good from bad ROIs.
3. Step through ROIs (arrow keys / Forward-Back), pressing **y**/**n** to label
   each; use **d** to flag doublets.
4. Adjust **CLim** / **FOV size** as needed for visibility.
5. Click **Export** to write the labeled table to CSV. Re-open later and use
   **Load** to continue from where you left off.

---

## Contents of this folder

- `annotateROI.mlapp` — the app.
- `uipickfiles.m` — file-selection dialog used by the app (its only in-repo
  dependency).
- `README.md` — this file.

Everything else the app uses (`tiffreadVolume`, `regionprops3`, `imdilate`,
`strel`, `prctile`, `uistyle`/`addStyle`, `tiledlayout`/`nexttile`, `listdlg`,
`readtable`/`writetable`, `contour`, …) ships with MATLAB and the two toolboxes
listed above, so no other project files are required.
