# SprayAR: Hand-Tracked Augmented Reality Viewer for Cold Spray Molecular Dynamics

<p align="center">
  <img src="https://img.shields.io/badge/WebXR-Immersive%20AR-blue?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Three.js-r160-black?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Meta%20Quest-Passthrough%20AR-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Hand%20Tracking-25%20Joints-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/OVITO-Data%20Pipeline-teal?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Source-LAMMPS%20Trajectories-yellow?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A pipeline that takes an existing <b>LAMMPS trajectory</b>, reads it through <b>OVITO</b>, and
  exports a self-contained <b>WebXR passthrough-AR viewer</b> for the Meta Quest -- no controllers,
  no native app, no marker cards. The exported point cloud floats in your real room at true
  human scale, colored by whichever per-atom field you choose (pressure, particle type, stress,
  coordination, and so on), and you interact with it entirely with your bare hands: pinch to
  grab, two-hand pinch to scale and rotate, or turn a palm face-up and let it rest there like a
  small holographic object sitting in your hand.
</p>

<img width="903" height="510" alt="ar cold" src="https://github.com/user-attachments/assets/0cc79cfb-255c-49a6-be40-69bd45c3921e" />

---

## Read This First: What This Pipeline Gives You, and What It Does Not

This is a **visualization and interaction layer**, not a solver. It does not run, validate, or
reinterpret your molecular dynamics simulation in any way -- it takes whatever positions and
per-atom fields your OVITO pipeline already computes and makes them walkable, graspable, and
palm-restable in a headset.

| This pipeline gives you | This pipeline does NOT give you |
|---|---|
| Direct, in-headset, hand-manipulable viewing of any existing LAMMPS/OVITO trajectory | New physics, a new simulation, or any solver of its own |
| Color mapping to any per-atom field OVITO can compute (pressure, stress, type, coordination, velocity, ...) | Automatic judgment of which field is scientifically meaningful -- you still choose `PROPERTY_NAME` |
| A real-world-scale, walk-around-able AR object anchored to your room | A permanently hosted viewer -- the output folder has to be re-hosted after every re-export |
| A lightweight binary export that keeps file size small regardless of atom count | Full-resolution playback -- frames are capped and subsampled for standalone headset hardware |
| Pinch-grab, two-hand scale/rotate, and palm-rest interaction, all built on the standard WebXR Hand Input API | A native Quest app, controller support, or any input method beyond bare-hand tracking |

If what you actually need is a rendering of *bonding physics*, *thermal history*, or anything
your OVITO pipeline itself doesn't already compute per atom, that has to be solved upstream --
this pipeline only visualizes whatever field you hand it.

---

## Why Binary Streaming Instead of Raw JSON

An early version of this pipeline embedded per-atom positions directly as JSON inside the HTML
file. That works for small trajectories, but text-encoded floating point numbers are roughly
three to five times larger than their raw binary equivalent, and cold-spray impact datasets can
easily carry hundreds of thousands of atoms per frame:

| Export format | Approx. size (200,000 atoms x 6 frames, 4 values/atom) |
|---|---|
| Raw JSON, embedded in the HTML | 60-90+ MB, in a single file the headset has to parse |
| Interleaved `[x, y, z, colorValue]` float32 binary, fetched per frame | ~18-19 MB, split across small files |

The exported HTML now fetches each frame's `.bin` file and reads it directly into a
`Float32Array` -- no JSON parsing on the critical path, and file size scales with atom count
instead of with atom count times text-encoding overhead.

---

## Why Auto-Detected Properties

Every LAMMPS dump names its per-atom quantities differently depending on what was computed
upstream -- `c_pe`, `v_pressure`, `Particle Type`, a custom compute, or nothing beyond position
at all. Hardcoding one property name ties the pipeline to a single dataset. Instead, the script
reads every particle property OVITO actually found in *your* file at runtime, prints the full
list to the console, and only falls back to `Particle Type` (always present in any LAMMPS dump)
if your requested field isn't there -- so the same script runs unmodified across a friction-stir
weld trajectory, a single-impact cold-spray run, and a post-impact deposition run, with only one
variable (`PROPERTY_NAME`) changing between them.

---

## Guaranteeing the Model Fits in Your Hand

Raw LAMMPS coordinates are not meters, and a naive 1:1 mapping would place the model anywhere
from invisible to the size of a building depending on your simulation's unit system. The
pipeline instead measures the trajectory's own bounding diagonal and rescales it to a fixed,
arm's-reach real-world size:

```
raw simulation extent (diagonal)     = allData.size            (arbitrary LAMMPS units)
desired real-world extent            = 0.35 m
initial scale factor                 = 0.35 / allData.size
```

Point size is compensated the same way. Three.js's `PointsMaterial.size` is not affected by an
object's own scale transform (a real quirk of how point-sprite sizing works in the renderer), so
without correction, shrinking the model down to real-world scale would leave the dots
disproportionately huge. The pipeline bakes the current scale into the material size manually
every time the model is rescaled, keeping dot size proportional whether you're looking at the
default view or one you've pinch-scaled up or down.

---

## Interaction Model

| Gesture | Effect |
|---|---|
| Single-hand pinch inside the model | Grab and reposition it |
| Two-hand pinch | Pull apart/together to scale; twist to rotate |
| Palm face-up beneath the model (Palm Mode) | Model hovers just above your palm and follows your hand |
| Frame slider / Play-Pause / First-Previous-Next-Last (2D panel) | Step through or auto-advance the exported frames |
| Recenter button | Resets position, scale, and rotation to the defaults |

Palm Mode and manual grab never fight each other: pinching always takes priority, and Palm Mode
resumes automatically the moment neither hand is pinching and a palm faces upward again.

---

## Reference Space: Why Height Means the Floor, Not Your Eyes

WebXR offers more than one coordinate origin convention. Requesting `local-floor` as a session
*feature* is not the same as telling three.js to actually use it -- without an explicit
`renderer.xr.setReferenceSpaceType('local-floor')` call, three.js defaults to the `local` space,
whose origin sits at roughly your eye level when the session starts. A model placed at `y = 1.2`
under that default would float 1.2 meters *above your eyes* rather than at chest height. This
pipeline sets the reference space explicitly, so a y-coordinate reliably means height above the
floor regardless of whether you're standing, sitting, or the headset started the session at an
odd angle.

---

## Repository Structure

```
sprayAR/
|
|-- export_ar_viewer.py            # Main pipeline script (edit INPUT_FILE / PROPERTY_NAME, run with OVITO)
|
|-- ar_output/                     # Generated by the script -- host this whole folder together
|   |-- <basename>_ar.html         # The WebXR AR viewer itself
|   |-- <basename>_frame_0.bin     # Interleaved [x, y, z, colorValue] float32 binary, one per exported frame
|   |-- <basename>_frame_1.bin
|   |-- ...
|
|-- README.md                      # This file
```

---

## How to Run

### Requirements

- OVITO Python environment: https://www.ovito.org
- Python 3.x: https://python.org
- A Meta Quest headset with Hand Tracking enabled (Settings > Movement Tracking)
- Any static file host reachable from the headset over HTTPS (or localhost via USB)

### Step 1 - Point the script at your trajectory

Open `export_ar_viewer.py` and set:

```python
INPUT_FILE = r"path\to\your_trajectory.lammpstrj"
PROPERTY_NAME = "v_pressure"
```

`PROPERTY_NAME` should match the exact per-atom field you want to color by. If it's not found,
the script falls back to `Particle Type` and prints every property it actually found in your
file, so you can copy the correct name in and re-run.

### Step 2 - Run the export

```bash
python export_ar_viewer.py
```

Console output confirms the detected properties, the frame indices chosen, the subsampling
applied (if any), and the default colormap picked for the field you're visualizing.

### Step 3 - Host the output folder

The exported HTML fetches its `.bin` files at runtime, so the whole `ar_output` folder has to be
hosted together, not just the HTML file on its own.

- Easiest: drag the entire `ar_output` folder onto a static-site drop host and use the link it
  gives you.
- For fast local iteration: serve the folder locally and reach it from the headset over USB
  (`adb reverse`) instead of re-uploading after every change.

### Step 4 - Open it on the Meta Quest

Open the built-in **Meta Quest Browser** (passthrough AR only works in Meta's own browser, not a
sideloaded Chrome) and navigate to your hosted link.

### Step 5 - Enter AR

Tap **Enter AR**, and accept both the passthrough camera and hand-tracking permission prompts.
The model appears roughly at chest height, about half a meter in front of wherever you were
facing when the session started, anchored to that spot in your room.

### Step 6 - Interact

Pinch to grab, two-hand pinch to scale/rotate, or turn a palm face-up beneath it to carry it
around. Use the 2D panel (aim and pinch like a laser pointer, or use the frame slider as a
larger, easier target) to step through frames or change the colormap.

---

## Common Errors and Fixes

| Symptom | Cause | Fix |
|---|---|---|
| Color looks the same everywhere | `PROPERTY_NAME` defaulted to a near-constant field such as `Particle Type` | Set `PROPERTY_NAME` to a continuous per-atom field (e.g. `v_pressure`) and re-run |
| Console prints `Property 'X' not found` | Exact column name/spelling/case differs in your dump file | Copy the exact name from the printed "Available particle properties" list |
| Model appears far above your head | Reference space defaulted to `local` (eye-level origin) instead of `local-floor` | Confirm `renderer.xr.setReferenceSpaceType('local-floor')` runs before `ARButton.createButton` |
| Palm Mode triggers when your palm faces down instead of up | `PALM_NORMAL_SIGN` convention is flipped for that hand | Flip the sign (`1` <-> `-1`) for that handedness |
| Blank page or a 404 fetching `.bin` files after hosting | The hosted URL points one folder level off from where the HTML actually landed | Check the host's file listing and adjust the URL to the exact path of `<basename>_ar.html` |
| Typing the hosted URL on the headset is painful | Auto-generated hosting subdomain is long and random | Rename the site to a short custom subdomain before opening it in the headset |
| Play/colormap buttons don't respond to a pinch | `dom-overlay` pinch-to-click behavior varies by headset/browser version | Aim more directly at the button and move closer, or use the frame slider (larger hit target) as a fallback |
| HTML/output folder balloons to hundreds of MB | `MAX_PARTICLES_PER_FRAME` or `NUM_FRAMES_TO_EXPORT` set too high for the dataset | Lower one or both and re-export |
| Playback looks like random flicker instead of continuous motion | Subsampling used a different random seed between runs | Keep `RNG_SEED` fixed across exports of the same dataset |

---

## Extending the Pipeline

| Extension | What to change |
|---|---|
| Different dataset or trajectory | `INPUT_FILE` |
| Different colored field | `PROPERTY_NAME` (auto-falls back and prints available options if the name is wrong) |
| More or fewer exported frames | `NUM_FRAMES_TO_EXPORT` |
| Performance vs. detail trade-off | `MAX_PARTICLES_PER_FRAME` |
| Bigger or smaller in-headset model | `DESIRED_WORLD_SIZE` |
| Where the model first appears | `DEFAULT_POSITION` |
| Palm Mode sensitivity or smoothness | `PALM_UP_THRESHOLD`, `PALM_FOLLOW_LERP` |
| Fixing inverted palm-up detection | `PALM_NORMAL_SIGN` |
| Default colormap for a given field | The keyword list checked against `PROPERTY_NAME` (defaults to `jet` for pressure/stress-like fields) |
| Reproducible subsampling across re-runs | `RNG_SEED` |

---

## Citation

If you use or adapt this pipeline, a citation template is provided below -- replace the
placeholders with your own author/DOI details once archived (e.g. via Zenodo):

```bibtex
@software{mishra_2026_sprayar,
  author    = {Mishra, A.},
  title     = {SprayAR: Hand-Tracked Augmented Reality Viewer for Cold Spray Molecular Dynamics},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21628146},
  url       = {https://doi.org/10.5281/zenodo.21628146}
}
```

> Mishra, A. (2026). *SprayAR: Hand-Tracked Augmented Reality Viewer for Cold Spray Molecular Dynamics* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21628146

---

## Author

**akshansh11**
GitHub: https://github.com/akshansh11

---

## License

<p align="center">
  <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
    <img alt="Creative Commons Licence" style="border-width:0; margin: 12px 0;"
      src="https://licensebuttons.net/l/by-nc/4.0/88x31.png"/>
  </a>
  <br/>
  <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
    Creative Commons Attribution-NonCommercial 4.0 International License
  </a>
</p>

This work is licensed under a **Creative Commons Attribution-NonCommercial 4.0 International License**.

You are free to:
- **Share** - copy and redistribute the material in any medium or format
- **Adapt** - remix, transform, and build upon the material

Under the following terms:
- **Attribution** - You must give appropriate credit to akshansh11 and link to this repository
- **NonCommercial** - You may not use the material for commercial purposes

Copyright (c) 2026 akshansh11. All rights reserved for commercial use.
