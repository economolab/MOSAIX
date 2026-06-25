# unmixV5

**Spectral unmixing of multichannel fluorescence images**

`unmixV5` is a MATLAB App Designer tool for cleaning up multichannel
fluorescence images in which signal from one fluorophore "bleeds" into a
spectrally adjacent channel.

## What this app does

- View the channels of a multi-page TIFF together as a color composite.
- Register (shift) channels relative to one another.
- Subtract spectral bleed-through from one channel into another
  (**spectral unmixing**).
- Optionally remove lipofuscin autofluorescence.
- Optionally apply a spatial (Gaussian) smoothing filter.
- Save/load settings and export the unmixed image stack, including batch
  processing of many files at once.

The right-hand image is a **live preview**. Almost every control updates this
preview so you can tune parameters by eye before exporting.

---

## Requirements

- **MATLAB** with **App Designer**, the **Image Processing Toolbox**, and the
  **Statistics and Machine Learning Toolbox**.
- The bundled helper files must be on the MATLAB path (keep them in the same
  folder as `unmixV5.mlapp`): `gaussFilt.m`, `normalize_local_contrast.m`,
  `qprctile.m`, `parseTiffChansPlanesInfo.m`, `uipickfiles.m`.

### NoRMCorre (optional — required only for non-rigid registration)

The **non-rigid registration** features — the **Apply non-rigid (normxcorr)**
button in the Image Shifts panel, and non-rigid correction during
**Export**/**Batch Process** — depend on the third-party
[**NoRMCorre**](https://github.com/flatironinstitute/NoRMCorre) toolbox
(`normcorre`, `NoRMCorreSetParms`, `apply_shifts`). It is **not** bundled here.

- If NoRMCorre is not on the MATLAB path, the first time you click
  **Apply non-rigid (normxcorr)** the app prompts you to select the NoRMCorre
  folder and adds it to the path for the session.
- All other features (channel display, image shifts, spectral unmixing,
  lipofuscin removal, spatial filter, export) work **without** NoRMCorre.

---

## 1. Getting started

### Launching & loading an image

- On start, the app opens a file picker. Select your image (`.tif`/`.tiff`,
  `.raw`, or `.prd`). The folder it lives in becomes the working folder.
- You are then prompted for two things:
  - **Enter plane to load** — which Z-plane / frame to display first
    (defaults to the middle plane).
  - **Image size** — leave as `inf` to load the full frame, or enter a smaller
    number to crop to a faster-to-preview region.
- All other TIFFs in the same folder are listed in the **File:** dropdown
  (top right) so you can switch between them later.

### Channel order (read this — it affects everything)

- The app orders channels by emission wavelength read from the file name. It
  looks for tokens ending in a 3-digit number between **250 and 999**
  (e.g. `...GFP510..._mCherry610...`). Channels are sorted low-to-high
  wavelength, and those tokens become the channel labels.
- If no wavelength numbers are found, channels are left in file order and
  labeled `Channel 1`, `Channel 2`, etc.
- This matters because the **Shift** and **Unmix** tables refer to channels by
  number (`1` = first/leftmost channel shown in the Channels panel).

---

## 2. Channels panel (display controls)

One column of controls appears per channel:

| Control | Purpose |
| --- | --- |
| **Name field** (top) | The channel label. Editable; purely cosmetic. |
| **Color swatch** | Click to open a color picker and set the channel's display color in the composite. |
| **Min / Max** | Black/white display limits (contrast). Lower values brighten dim signal. Editable only when **Auto** is OFF. |
| **Auto** | When ON, Min/Max are set automatically from image percentiles (~0.1%–99.995%) and locked. Turn OFF to type your own limits. |
| **On** | Show/hide this channel in the composite. First three channels are ON by default. Affects preview only, not exported data. |

- **HiLo LUT** (checkbox, right of the channels): overlays blue on every pixel
  that is exactly zero in all displayed channels — useful for spotting
  saturated/clipped or empty regions while setting limits.

> **Note:** The top **Plane:** field and the **Auto All** button are
> labels/placeholders in this build — the active plane is chosen via the
> dialog when you load or switch files, and **Auto All** has no action wired to
> it. Use the per-channel **Auto** checkboxes instead.

---

## 3. Image Shifts panel — registering channels to each other

Use this when channels are spatially offset (chromatic shift, camera
misalignment). Bleed-through subtraction only works well once channels are
aligned.

**Shift table columns:**

| Column | Meaning |
| --- | --- |
| **Channel** | The channel number to move. |
| **X shift** | Horizontal shift (positive = right). |
| **Y shift** | Vertical shift (positive = down). |

Add one row per channel you want to move; leave the others at `0`.

- Whole-number shifts move the channel by that many pixels.
- If **any** shift value is non-integer (e.g. `1.5`), the app switches to a
  sub-pixel Fourier (FFT) shift for all channels instead.

**Buttons:**

- **Apply non-rigid (normxcorr):** runs NoRMCorre non-rigid registration to
  warp channels 1 and 2 onto each other (for local, non-uniform distortion).
  If the NoRMCorre toolbox isn't on the MATLAB path, you'll be asked to point
  to its folder once. This modifies the loaded stack in memory.
- **Revert non-rigid:** reloads the current file from disk, discarding the
  non-rigid warp (table shifts and unmix settings are re-applied on reload).

After editing the shift table, click **Update image** (section 6) to see the
result.

---

## 4. Spectral Unmixing panel — the core bleed-through subtraction

Each row defines one subtraction. Add as many rows as you need.

| Column | Meaning |
| --- | --- |
| **Unmix from** | The SOURCE channel (the one bleeding/contaminating). |
| **Unmix to** | The TARGET channel you want to clean up. |
| **Offset** | Background level subtracted from the source before use. |
| **Scale** | The bleed-through coefficient (fraction of source to subtract). |
| **Enable** | Tick to apply this row; untick to ignore it. |

The math applied to each enabled row is:

```
cleaned = target - max( (source - Offset) * Scale , 0 )
```

i.e. the source image is background-subtracted (`Offset`), negatives clamped
to zero, scaled by `Scale`, then subtracted from the target.

### How to tune by hand

1. Pick a region where you **know** only the source fluorophore is present (so
   any signal in the target there is pure bleed-through).
2. Set **Offset** to roughly that channel's background value.
3. Increase **Scale** until the bleed-through in the target disappears in the
   preview, without creating dark "holes" (over-subtraction).
4. Tick **Enable** and click **Update image**.

### Auto Estimate (button)

Automatically fills the table and estimates **Scale** for each adjacent
channel pair using linear regression between channels (both directions:
channel `i` into `i+1`, and `i` into `i-1`). Rows whose estimate is
essentially zero are left disabled.

> **Important:** for a good estimate it needs a background region defined.
> First go to the **Lipo options** panel, tick **Show ROIs**, and drag the
> background rectangle over a true-background area (see section 5).

**Order of operations:** when the image is rebuilt, processing runs as
`refresh → lipofuscin removal → channel shifts → spatial filter → spectral
unmixing`. So shifts and lipo removal happen before unmixing.

---

## 5. Lipo options panel — lipofuscin removal & background ROI

Lipofuscin is broadband autofluorescence that appears in many channels at
once. This panel models it and subtracts it.

- **Show ROIs** (checkbox): shows two tools on the image:
  - a **rectangle** (background ROI) — drag/resize over a region of pure
    background (no real signal). Also used by **Auto Estimate**.
  - a **crosshair** — place it directly on a bright lipofuscin spot.
- The small plot shows the modeled background spectrum and the lipofuscin
  amplitude spectrum across channels.
- **Channels to build model from:** channel numbers used to detect where
  lipofuscin is (e.g. `1 2 3 4 5`).
- **Channels to remove lipo from:** channel numbers to subtract it out of.
- **Update lipo model** (button): recomputes the lipofuscin background and
  amplitude from the current crosshair (a lipo spot) and rectangle
  (background) positions. Do this whenever you move the ROIs.
- **Remove lipofuscin** (checkbox): when ticked, the modeled lipofuscin is
  subtracted as part of the processing pipeline.

**Typical use:** tick **Show ROIs**, position the rectangle on background and
the crosshair on a lipofuscin granule, set the model/remove channel lists,
click **Update lipo model**, tick **Remove lipofuscin**, then **Update image**.

### Spatial filter (bottom of this panel)

- **Spatial Filter** (checkbox): apply a Gaussian smoothing filter to reduce
  noise before/around unmixing.
- **width (pixels):** the Gaussian filter radius (default `1.25`). Larger =
  more smoothing.

---

## 6. Applying, saving, and exporting (right-hand buttons)

| Button | Action |
| --- | --- |
| **Update image** | Re-runs the whole pipeline (lipo removal, shifts, spatial filter, unmixing) on the current plane and refreshes the preview. Click after changing any table or option. |
| **Save Settings** | Saves all current settings (labels, colors, limits, enable/auto states, shift table, unmix table, lipo model, spatial filter) to a `.mat` file. |
| **Load Settings** | Loads a previously saved `.mat` settings file and re-applies it to the current image. Use this to apply the same recipe to new images. |
| **Export Image** | Processes **every plane** of the current file with current settings and writes a multi-page 16-bit TIFF named `<originalname>-unmixed.tif` in the source folder. Channels are written back in original (file) order. (If non-rigid registration was applied, it is re-run per plane during export.) |
| **Batch Process** | Applies current settings to every file selected in the **File:** dropdown and exports an `-unmixed.tif` for each. `Done Batch Processing!!!!` prints to the MATLAB console when finished. |
| **Composite Fig** | Opens a separate MATLAB figure showing the merged composite plus one panel per channel (each in its assigned color) — handy for figures/QC. |

---

## 7. Suggested workflow

1. Launch and select your image; choose a plane and full size (`inf`).
2. In the **Channels** panel, confirm channel order/labels, set display
   colors, and adjust contrast (use **Auto**, or turn it off for manual
   Min/Max).
3. If channels are misaligned, fill the **Image Shifts** table and
   **Update image**.
4. Open **Lipo options**, tick **Show ROIs**, place the background rectangle
   and the lipofuscin crosshair, set the channel lists, and
   **Update lipo model**. Tick **Remove lipofuscin** if needed.
5. In **Spectral Unmixing**, either click **Auto Estimate** (with the
   background ROI placed) or add rows manually; tune Offset/Scale and Enable
   each row.
6. Click **Update image** and inspect the preview (toggle channels On/Off and
   use **HiLo LUT** to check for over-subtraction).
7. **Save Settings** so you can reuse the recipe.
8. **Export Image** for this file, or **Batch Process** to apply to all files
   in the dropdown.

---

## Tips / troubleshooting

- **Over-subtraction** shows up as dark holes where real signal used to be.
  Lower the **Scale** (or raise the **Offset**) for that unmix row.
- Channels must be **aligned before unmixing** for the subtraction to be clean.
- **Auto Estimate** and good lipofuscin removal both depend on a correctly
  placed background ROI — take the time to position it on true background.
- Tables reference channels by **number** (`1` = first channel in the Channels
  panel), which follows the wavelength-sorted order described in section 1.
- Changes don't take effect until you click **Update image**.
