# SprayAR: Hand-Tracked Augmented Reality Viewer for Cold Spray Molecular Dynamics

<p align="center">
  <img src="https://img.shields.io/badge/WebXR-Immersive%20AR-blue?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Three.js-r160-black?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Meta%20Quest-Passthrough%20AR-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Hand%20Tracking-25%20Joints-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Rendering-Sphere%20Impostors-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/OVITO-Data%20Pipeline-teal?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Source-LAMMPS%20Trajectories-yellow?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A pipeline that takes an existing <b>LAMMPS trajectory</b>, reads it through <b>OVITO</b>, and
  exports a self-contained <b>WebXR passthrough-AR viewer</b> for the Meta Quest -- no controllers,
  no native app, no marker cards. Atoms are drawn as <b>shaded, depth-correct spheres</b> with
  ambient occlusion, close to the look of an OVITO Tachyon render, and the full impact sequence
  plays back inside AR. The model floats in your real room at arm's-reach scale, colored by
  whichever per-atom field you choose (pressure, particle type, stress, coordination, and so on),
  and you interact with it entirely with your bare hands: pinch to grab, two-hand pinch to scale
  and rotate, turn a palm face-up to let it rest in your hand, and poke a floating button panel
  to play, pause, and step through frames.
</p>

<img width="903" height="510" alt="ar cold" src="https://github.com/user-attachments/assets/0cc79cfb-255c-49a6-be40-69bd45c3921e" />

---

## What's New in v2

| Area | v1 | v2 |
|---|---|---|
| Atom rendering | Flat, unlit `PointsMaterial` squares | GPU sphere impostors: per-pixel normals, diffuse + specular lighting, darkened silhouettes, correct per-pixel depth so spheres intersect properly |
| Shading depth cues | None | Precomputed ambient occlusion: atoms in crevices and folds are darker |
| Which atoms are exported | Random subsample of the whole volume (speckled, interior atoms showing through) | Surface atoms only, found by coordination analysis; every visible atom is kept |
| Color scaling | Absolute min/max per frame (outliers squash everything into one color) | Percentile-clipped range, fixed across all frames, adjustable live in the viewer |
| Colormap switching | Re-fetched and re-colored every frame | Instant, via a GPU lookup texture |
| Atom size | Arbitrary fraction of the box size | Real atomic radius in Angstrom, per type, with a live size slider |
| Frame selection | Fixed number of evenly spaced frames | Start / end / stride, or evenly spaced; several dump files can be joined into one sequence |
| Model orientation | Simulation axes used as-is (LAMMPS Z ended up horizontal) | Configurable up axis so the particle is above the substrate; switchable live |
| Playback inside AR | Not possible (2D buttons are not visible in Quest AR; timer-based playback stalled) | Floating in-AR button panel you press with your fingertip; playback driven by the render loop |
| Memory on the headset | All frames loaded at once | Limited frame cache with look-ahead prefetching, so long animations fit in headset memory |
| Frame labels | Frame index | Real LAMMPS timestep |
| Diagnostics | None | `check_frames.py` shows whether a trajectory actually contains the approach and impact |

---

## Read This First: What This Pipeline Gives You, and What It Does Not

This is a **visualization and interaction layer**, not a solver. It does not run, validate, or
reinterpret your molecular dynamics simulation in any way -- it takes whatever positions and
per-atom fields your OVITO pipeline already computes and makes them walkable, graspable, and
palm-restable in a headset.

| This pipeline gives you | This pipeline does NOT give you |
|---|---|
| Direct, in-headset, hand-manipulable viewing of any existing LAMMPS/OVITO trajectory | New physics, a new simulation, or any solver of its own |
| Shaded sphere rendering with ambient occlusion, close to an OVITO Tachyon image | Ray-traced shadows, true global illumination, or pixel-identical OVITO renders |
| Color mapping to any per-atom field OVITO can compute (pressure, stress, type, coordination, velocity, ...) | Automatic judgment of which field is scientifically meaningful -- you still choose `PROPERTY_NAME` |
| Animated playback of the impact sequence, controllable inside AR | Frames your LAMMPS run never dumped -- the trajectory must contain the approach and impact |
| A real-world-scale, walk-around-able AR object anchored to your room | A permanently hosted viewer -- the output folder has to be re-hosted after every re-export |
| A compact binary export of surface atoms only | Views of the interior -- interior atoms are removed (set `SURFACE_ONLY = False` to keep them) |
| Pinch-grab, two-hand scale/rotate, palm-rest, and fingertip buttons, all on the standard WebXR Hand Input API | A native Quest app, controller support, or any input method beyond bare-hand tracking |

---

## How the Rendering Works

### Sphere impostors

Real sphere meshes are far too expensive for hundreds of thousands of atoms on standalone headset
hardware. Each atom is instead a single GPU point sprite, sized in the vertex shader to the
atom's projected diameter on screen:

```
pointSize_px = viewportHeight_px * P[1][1] * radius_view / depth
```

The fragment shader then turns each square sprite into a sphere:

1. Pixels outside the unit circle are discarded.
2. A surface normal is reconstructed from the pixel's position inside the circle.
3. The sphere surface point is projected back to write a **correct depth** for that pixel, so
   neighboring atoms intersect like real spheres instead of overlapping as flat discs.
4. Diffuse and specular lighting from a head-mounted light, a darkened silhouette (similar to
   Tachyon's edge look), and the atom's ambient-occlusion factor are applied.

Because the radius is computed from the model's world transform, atoms stay correctly sized when
you pinch-scale the model -- no manual point-size correction is needed.

### Surface-only export

Interior atoms can never be seen, so exporting them only costs file size and headset performance.
OVITO's `CoordinationAnalysisModifier` counts neighbors within `COORD_CUTOFF`. Atoms with fewer
than `BULK_COORDINATION` neighbors are kept as surface atoms. Periodic boundaries are switched off
for this step so the side faces of the substrate slab count as surfaces.

| Lattice | `COORD_CUTOFF` | `BULK_COORDINATION` |
|---|---|---|
| FCC (Cu, Al, Ni, ...) | Between 1st and 2nd neighbor shell (Cu: ~3.0 A) | 10 (12 is a perfect lattice; 10 tolerates impact heating and deformation) |
| BCC (Fe, ...) | Just above the 2nd shell (Fe: ~3.2 A) | 12 (14 is a perfect lattice) |

### Ambient occlusion

For each exported atom, the number of atoms within `AO_RADIUS` is counted with a k-d tree
(`scipy`). Atoms in concave regions -- folds, jet roots, the gap where the particle meets the
substrate -- have more neighbors and are darkened. This is computed once at export time and costs
nothing in the headset.

### Color range and colormaps

The color range is taken from percentiles (default 2nd-98th) of the property across all exported
frames, so a handful of extreme atoms cannot squash the rest into a single color, and the same
color means the same value in every frame. Raw property values are stored in the binary files and
normalized on the GPU, so the range can be changed live in the viewer, and colormaps are switched
instantly through a 256-entry lookup texture.

---

## Binary Frame Format

Each frame is one `.bin` file of interleaved `float32` values, six per atom:

```
[x, y, z, value, radius, ao]   x, y, z centered on the first frame's centroid
```

The viewer reads each file directly into a `Float32Array` and hands it to the GPU without
copying. File size is about **24 bytes per exported atom**:

```
frame size  ~ 24 B x surface atoms
total size  ~ 24 B x surface atoms x number of frames
```

Surface extraction typically removes most of the atoms of a bulk substrate, which keeps frames
small even for million-atom systems. Only the frames around the current one are held in headset
memory (`MAX_CACHED_FRAMES`), and the next few are loaded ahead (`PREFETCH_AHEAD`), so the total
number of frames is limited by hosting and download time rather than headset RAM.

---

## Why Auto-Detected Properties

Every LAMMPS dump names its per-atom quantities differently depending on what was computed
upstream -- `c_pe`, `v_pressure`, `Particle Type`, a custom compute, or nothing beyond position
at all. The script reads every particle property OVITO actually found in *your* file at runtime,
prints the full list to the console, and only falls back to `Particle Type` (always present in any
LAMMPS dump) if your requested field isn't there. The same script runs unmodified across a
friction-stir weld trajectory, a single-impact cold-spray run, and a post-impact deposition run,
with only `PROPERTY_NAME` changing between them. The default colormap is `jet` for pressure- and
stress-like fields and `viridis` otherwise.

---

## Guaranteeing the Model Fits in Your Hand

Raw LAMMPS coordinates are not meters. The pipeline measures the bounding-box diagonal of the
exported atoms **across all exported frames** -- including the sprayed particle at its starting
height -- and rescales it to a fixed arm's-reach size:

```
raw simulation extent (diagonal)     = size over all frames       (Angstrom)
desired real-world extent            = 0.35 m                     (DESIRED_WORLD_SIZE)
initial scale factor                 = 0.35 / size
```

The model therefore stays in view during the whole approach, impact, and deformation.

---

## Orientation: Particle on Top, Substrate Below

LAMMPS simulations usually use **Z** as the vertical (impact) axis, while WebXR and three.js use
**Y** as "up". Without correction, the substrate appears standing on its edge. The viewer places
the data inside an orientation group that maps the simulation's up axis onto the world's up axis:

| `UP_AXIS` | Use when |
|---|---|
| `+z` (default) | The particle starts at high Z and moves toward lower Z |
| `-z` | The particle starts at low Z (impact from below) |
| `+x`, `-x`, `+y`, `-y` | The impact direction runs along another axis |

Use the **Up axis** button (desktop panel or in-AR panel) to cycle through the options until the
particle is above the substrate, then put that value into `UP_AXIS` so it becomes the default.
Hand rotation and Recenter act on the outer model group, so they never undo the orientation.

---

## Interaction Model

### Hand gestures (in AR)

| Gesture | Effect |
|---|---|
| Single-hand pinch inside the model | Grab and reposition it |
| Two-hand pinch | Pull apart/together to scale; twist to rotate |
| Palm face-up beneath the model (Palm Mode) | Model hovers just above your palm and follows your hand |

Palm Mode and manual grab never fight each other: pinching always takes priority, and Palm Mode
resumes automatically the moment neither hand is pinching and a palm faces upward again.

### In-AR button panel

The 2D page controls are not visible inside Quest passthrough AR, so the viewer adds a floating
3D panel at chest height, about 40 cm in front of where you are facing when AR starts. Press a
button by touching it with your **index fingertip**:

| Button | Effect |
|---|---|
| `<<` / `>>` | Previous / next frame |
| `Play` / `Pause` | Start or stop the animation |
| `Up` | Cycle the up axis |
| `Reset` | Return the model to its default position, scale, and rotation |

A label above the buttons shows the current frame and LAMMPS timestep. If you walk more than
1.5 m away, the panel moves back in front of you.

### Desktop preview controls

| Control | Effect |
|---|---|
| Drag / scroll | Rotate / zoom |
| Arrow keys, Space | Step frames, play/pause |
| First / Previous / Play / Next / Last / slider | Frame navigation |
| Colormap dropdown, Max/Min fields, Auto range | Color settings (instant) |
| Atom size slider | Scale all atom radii |
| Up axis | Change the vertical axis |
| Recenter, Palm Mode | Reset the model; toggle palm-rest |

Settings chosen on the page before tapping **Enter AR** (colormap, range, atom size) carry over
into the AR session.

---

## Reference Space: Why Height Means the Floor, Not Your Eyes

WebXR offers more than one coordinate origin convention. Requesting `local-floor` as a session
*feature* is not the same as telling three.js to actually use it -- without an explicit
`renderer.xr.setReferenceSpaceType('local-floor')` call, three.js defaults to the `local` space,
whose origin sits at roughly your eye level when the session starts. A model placed at `y = 1.2`
under that default would float 1.2 meters *above your eyes* rather than at chest height. This
pipeline sets the reference space explicitly, so a y-coordinate reliably means height above the
floor. Foveated rendering is also switched off (`setFoveation(0)`) so atoms at the edge of your
view stay sharp.

---

## Repository Structure

```
sprayAR/
|
|-- export_ar.py                   # Main pipeline script (edit SETTINGS, run with OVITO's Python)
|-- check_frames.py                # Diagnostic: does the trajectory contain the approach and impact?
|
|-- ar_output/                     # Generated by the script -- host this whole folder together
|   |-- <basename>_ar.html         # The WebXR AR viewer itself
|   |-- <basename>_frame_0.bin     # Interleaved [x, y, z, value, radius, ao] float32, one per frame
|   |-- <basename>_frame_1.bin
|   |-- ...
|
|-- README.md                      # This file
```

---

## How to Run

### Requirements

- Python 3.11 (e.g. via Anaconda): https://www.anaconda.com
- OVITO Python module (`pip install ovito`): https://www.ovito.org
- `numpy`, and `scipy` for ambient occlusion (optional but recommended)
- A Meta Quest headset with Hand Tracking enabled (Settings > Movement Tracking)
- Any static file host reachable from the headset over HTTPS (or localhost via USB)

### Step 0 - One-time environment setup (Anaconda Prompt)

```bash
conda create -n ovito_env python=3.11 -y
conda activate ovito_env
pip install ovito numpy scipy
python -c "import ovito; print(ovito.version_string)"
```

Run `conda create` only once; running it again replaces the environment and removes OVITO. On an
unstable connection, use `pip install --no-cache-dir --timeout 120 --retries 10 ovito numpy scipy`,
or download the matching `.whl` from https://pypi.org/project/ovito/#files and install it locally.

In every new Anaconda Prompt session afterwards, only run `conda activate ovito_env`.

### Step 1 - Check what the trajectory contains

Set `INPUT_FILE` in `check_frames.py` and run:

```bash
python check_frames.py
```

It prints the frame count, timestep, atom count, and maximum Z for sampled frames. If max Z drops
over the frames, the particle approach is in the file. If it is constant from the first frame,
the file only holds the deposited state, and you need the dump from the impact run.

### Step 2 - Configure the export

Open `export_ar.py` and edit the SETTINGS block. The most important options:

```python
INPUT_FILE = r"path\to\your_trajectory.lammpstrj"
MULTIPLE_FRAMES = False        # True when INPUT_FILE is a wildcard like deposition_*.lammpstrj
PROPERTY_NAME = "v_pressure"

FRAME_START = 0                # export window
FRAME_END = None
FRAME_STRIDE = 1
NUM_FRAMES_TO_EXPORT = None    # or a number for evenly spaced frames

TYPES_TO_SHOW = None           # e.g. [2] to show only the sprayed particle
UP_AXIS = '+z'                 # simulation axis that points up in AR

COORD_CUTOFF = 3.0             # surface extraction, see the lattice table above
BULK_COORDINATION = 10
DEFAULT_RADIUS = 1.28          # Angstrom, ~half the nearest-neighbor distance
```

### Step 3 - Run the export

```bash
D:
cd "D:\path\to\your\project"
python export_ar.py
```

For each frame the console prints the total atoms, the surface atoms found, and the atoms
exported, followed by the automatic color range and the total output size. A
`WARNING ... subsampled` line means surface extraction kept too many atoms; lower
`BULK_COORDINATION` by one or raise `COORD_CUTOFF` slightly.

### Step 4 - Preview on your PC

```bash
cd ar_output
python -m http.server 8000
```

Open `http://localhost:8000/<basename>_ar.html` (default: `deposition_ar.html`) and press
Ctrl+F5. Do not double-click the HTML file: browsers block loading the `.bin` files from
`file://`. The server window should show `200` responses for the `.bin` files.

### Step 5 - Host the output folder

The HTML fetches its `.bin` files at runtime, so the whole `ar_output` folder has to be hosted
together.

- **Netlify Drop (easiest):** log in at https://app.netlify.com/drop and drag the whole
  `ar_output` folder onto the page. Rename the HTML to `index.html` first so the viewer opens at
  the site root, and give the site a short name under *Site configuration > Change site name*.
  After a re-export, drag the new folder onto the *Deploys* tab of the same site.
- **USB, no upload:** keep the local server running, enable developer mode on the Quest, run
  `adb reverse tcp:8000 tcp:8000`, and open `http://localhost:8000/<basename>_ar.html` in the
  Quest Browser.

Anyone with a hosted link can open it; use the USB route for unpublished data.

### Step 6 - Open it on the Meta Quest

Open the built-in **Meta Quest Browser** (passthrough AR only works in Meta's own browser) and
navigate to your link. Choose the colormap, color range, and atom size on the page.

### Step 7 - Enter AR

Tap **Enter AR** and accept the passthrough and hand-tracking permission prompts. The model
appears at chest height, about half a meter in front of you, and the button panel appears below
your line of sight.

### Step 8 - Interact

Pinch to grab, two-hand pinch to scale/rotate, or turn a palm face-up to carry the model. Press
`Play` on the floating panel with your fingertip to watch the impact.

---

## Common Errors and Fixes

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: No module named 'ovito'` | Wrong environment active, or the environment was recreated | `conda activate ovito_env`; if it was recreated, run `pip install ovito numpy scipy` again |
| pip download stops with "Connection interrupted" / `WinError 32` | Unstable network; antivirus locks the partial file | Use `--no-cache-dir --timeout 120 --retries 10`, or install a browser-downloaded `.whl` |
| "Failed to fetch" in the viewer | HTML opened by double-click (`file://`) | Serve the folder and open it via `http://localhost:8000/...` |
| Server shows `404` for `.bin` files | Server started in the wrong folder, or an old tab is still requesting files | Start the server inside `ar_output`; close old tabs and reload with Ctrl+F5 |
| `cd ar_output` fails | The export did not run or failed | Run `export_ar.py` and wait for "AR VIEWER CREATED" |
| Only the final deposited state is shown | The trajectory has no approach frames, or too few frames were exported | Run `check_frames.py`; export a denser frame window; join impact and post-impact dumps with a wildcard and `MULTIPLE_FRAMES = True` |
| Substrate stands on its edge / particle beside it | Simulation up axis differs from the viewer's Y-up | Press **Up axis** until the particle is on top, then set `UP_AXIS` |
| Animation does not play inside AR | 2D page buttons are not visible in Quest AR | Use the floating panel's `Play` button with your fingertip |
| Floating panel buttons do not react | Fingertip not touching the button plane, or hand tracking lost | Touch the button face directly with the index fingertip; keep your hand in view of the headset |
| `WARNING ... subsampled` on many frames, speckled look | Most atoms counted as surface (heated/deformed lattice) | Lower `BULK_COORDINATION` (e.g. 10 -> 9) or raise `COORD_CUTOFF` slightly |
| Holes in the surface | Surface criterion too strict | Raise `BULK_COORDINATION` by one or lower `COORD_CUTOFF` |
| Box sides missing / hollow-looking substrate | Periodic boundaries counted the side faces as interior | Keep the built-in PBC switch-off (`disable_pbc`) in the pipeline |
| Almost all atoms one color | Range dominated by outliers, or a near-constant field | Use the Max/Min fields or `COLOR_RANGE`; try `SYMMETRIC_RANGE = True` for pressure; check `PROPERTY_NAME` |
| Console prints `Property 'X' not found` | Exact column name/spelling/case differs in your dump file | Copy the exact name from the printed "Available particle properties" list |
| Gaps between atoms, or atoms overlapping heavily | Radius does not match the lattice | Use the Atom size slider, then set `DEFAULT_RADIUS` or `TYPE_RADII` |
| Model appears far above your head | Reference space defaulted to `local` instead of `local-floor` | Confirm `renderer.xr.setReferenceSpaceType('local-floor')` runs before `ARButton.createButton` |
| Palm Mode triggers when your palm faces down | `PALM_NORMAL_SIGN` convention is flipped for that hand | Flip the sign (`1` <-> `-1`) for that handedness |
| Stuttering in the headset | Too many atoms per frame | Narrow the frame window, lower `MAX_PARTICLES_PER_FRAME`, or tighten surface extraction |
| Output folder very large | Many frames x many surface atoms | Export only the impact window densely, use a larger `FRAME_STRIDE` elsewhere |

---

## Extending the Pipeline

| Extension | What to change |
|---|---|
| Different dataset or trajectory | `INPUT_FILE` (wildcard + `MULTIPLE_FRAMES = True` for several files) |
| Different colored field | `PROPERTY_NAME` |
| Which frames are exported | `FRAME_START`, `FRAME_END`, `FRAME_STRIDE`, `NUM_FRAMES_TO_EXPORT` |
| Show only the sprayed particle or the substrate | `TYPES_TO_SHOW` |
| Vertical axis in AR | `UP_AXIS` |
| Surface detection for another lattice | `COORD_CUTOFF`, `BULK_COORDINATION` |
| Keep interior atoms (e.g. for future slicing) | `SURFACE_ONLY = False` |
| Atom radii | `DEFAULT_RADIUS`, `TYPE_RADII` |
| Strength of crevice shading | `AO_ENABLED`, `AO_RADIUS`, `AO_STRENGTH` |
| Fixed or symmetric color range | `COLOR_RANGE`, `CLIP_PERCENTILES`, `SYMMETRIC_RANGE` |
| Playback speed | `PLAY_INTERVAL_MS` |
| Headset memory vs. smoothness | `MAX_CACHED_FRAMES`, `PREFETCH_AHEAD` |
| Performance safety cap | `MAX_PARTICLES_PER_FRAME` |
| Bigger or smaller in-headset model | `DESIRED_WORLD_SIZE` (in the HTML template) |
| Where the model first appears | `DEFAULT_POSITION` (in the HTML template) |
| Lighting direction | `uLightDir` uniform (in the HTML template) |
| Palm Mode sensitivity or smoothness | `PALM_UP_THRESHOLD`, `PALM_FOLLOW_LERP` |
| Fixing inverted palm-up detection | `PALM_NORMAL_SIGN` |
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
