# -*- coding: utf-8 -*-
"""
Cold spray deposition -> WebXR AR viewer with Tachyon-like sphere rendering

Run from the Anaconda Prompt:
    conda activate ovito_env
    D:
    cd "D:\\Atomisitc modeling research\\Electrostatic cold spray deposition"
    python export_ar.py

Then preview:
    cd ar_output
    python -m http.server 8000
    -> open http://localhost:8000/<BASENAME>_ar.html
"""
import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI.*')
import ovito._extensions.pyscript

from ovito.io import import_file
from ovito.modifiers import (CoordinationAnalysisModifier,
                             ExpressionSelectionModifier,
                             DeleteSelectedModifier)
import numpy as np
import os
import json

try:
    from scipy.spatial import cKDTree
except ImportError:
    cKDTree = None

print("=" * 70)
print("COLD SPRAY DEPOSITION - AR VIEWER EXPORT (sphere impostors)")
print("=" * 70)

# ----------------------------------------------------------------------------
# SETTINGS
# ----------------------------------------------------------------------------
INPUT_FILE = r"D:\Atomisitc modeling research\Electrostatic cold spray deposition\deposition_postimpact.lammpstrj"
# To combine several dump files into one sequence (e.g. impact + post-impact),
# use a wildcard and set MULTIPLE_FRAMES = True. Files must sort in time order by name:
#   INPUT_FILE = r"D:\...\deposition_*.lammpstrj"
MULTIPLE_FRAMES = False

OUTPUT_DIR = os.path.join(os.path.dirname(INPUT_FILE), "ar_output")
BASENAME = "deposition"
TITLE = "Cold Spray Deposition"

PROPERTY_NAME = "v_pressure"

# Which frames to export
FRAME_START = 0              # first frame index
FRAME_END = None             # last frame index (None = last frame in file)
FRAME_STRIDE = 1             # export every Nth frame between START and END
NUM_FRAMES_TO_EXPORT = None  # or a number (e.g. 60) for evenly spaced frames instead of STRIDE

# Which particle types to show. None = all. E.g. [2] to show only the sprayed
# particle and hide the substrate.
TYPES_TO_SHOW = None

# Only export surface atoms (interior atoms are invisible anyway).
#   FCC (Cu, Al, Ni...): cutoff between 1st and 2nd shell (Cu: ~3.0 A), bulk = 12
#   BCC (Fe...):        cutoff just above 2nd shell (Fe: ~3.2 A),      bulk = 14
SURFACE_ONLY = True
COORD_CUTOFF = 3.0
BULK_COORDINATION = 10     # 12 = perfect FCC; 10 tolerates thermal/deformation disorder

# Atom radius in Angstrom (~ half the nearest-neighbour distance looks like OVITO)
DEFAULT_RADIUS = 1.28
TYPE_RADII = {}            # per-type override, e.g. {1: 1.28, 2: 1.43}

# Ambient occlusion (needs scipy). Atoms in crevices get darker.
AO_ENABLED = True
AO_RADIUS = 7.0            # Angstrom, ~3x nearest-neighbour distance
AO_STRENGTH = 0.55         # 0 = off, 1 = very strong

# Color range: None = automatic percentile clipping over all exported frames
COLOR_RANGE = None         # e.g. (-50000, 50000)
CLIP_PERCENTILES = (2, 98)
SYMMETRIC_RANGE = False    # True -> range centred on zero (nice for pressure)

MAX_PARTICLES_PER_FRAME = 400_000   # safety cap after surface extraction
RNG_SEED = 42

# Which simulation axis points 'up' in AR (particle on top, substrate below).
# One of '+z', '-z', '+x', '-x', '+y', '-y'. Can also be changed live with the 'Up axis' button.
UP_AXIS = '+z'

# Viewer playback
PLAY_INTERVAL_MS = 100     # time per frame during playback (100 ms = 10 fps)
MAX_CACHED_FRAMES = 12     # frames kept in headset memory at once
PREFETCH_AHEAD = 3         # frames loaded ahead of the current one
# ----------------------------------------------------------------------------

if MULTIPLE_FRAMES:
    pipeline = import_file(INPUT_FILE, multiple_frames=True)
else:
    pipeline = import_file(INPUT_FILE)

if TYPES_TO_SHOW:
    expr = " && ".join(f"ParticleType != {t}" for t in TYPES_TO_SHOW)
    pipeline.modifiers.append(ExpressionSelectionModifier(expression=expr))
    pipeline.modifiers.append(DeleteSelectedModifier())

if SURFACE_ONLY:
    # Turn off PBC so that the side faces of the slab count as surfaces
    def disable_pbc(frame, data):
        data.cell_.pbc = (False, False, False)
    pipeline.modifiers.append(disable_pbc)
    pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff=COORD_CUTOFF, number_of_bins=10))

total_frames = pipeline.source.num_frames
print(f"Total frames: {total_frames}")

data0 = pipeline.compute(0)
available_props = list(data0.particles.keys())
print(f"Available particle properties: {available_props}")
if PROPERTY_NAME not in available_props:
    fallback = "Particle Type" if "Particle Type" in available_props else available_props[0]
    print(f"Property '{PROPERTY_NAME}' not found -> using '{fallback}' instead.")
    PROPERTY_NAME = fallback

if SURFACE_ONLY:
    print(f"Surface extraction ON (cutoff {COORD_CUTOFF} A, bulk coordination {BULK_COORDINATION})")
if AO_ENABLED and cKDTree is None:
    print("scipy not installed -> ambient occlusion disabled (pip install scipy to enable)")


def get_radii(data, n):
    radii = np.full(n, DEFAULT_RADIUS, dtype=np.float32)
    try:
        ptypes = data.particles.particle_types
        ids = np.asarray(ptypes)
        for t in ptypes.types:
            if t.id in TYPE_RADII:
                radii[ids == t.id] = TYPE_RADII[t.id]
    except Exception:
        pass
    return radii


def compute_ao(all_pos, query_pos):
    if not AO_ENABLED or cKDTree is None or len(query_pos) == 0:
        return np.ones(len(query_pos), dtype=np.float32)
    tree = cKDTree(all_pos)
    counts = tree.query_ball_point(query_pos, r=AO_RADIUS, return_length=True, workers=-1).astype(np.float64)
    lo, hi = np.percentile(counts, [5, 95])
    occl = np.clip((counts - lo) / (hi - lo), 0, 1) if hi > lo else np.zeros_like(counts)
    return (1.0 - AO_STRENGTH * occl).astype(np.float32)


end = total_frames - 1 if FRAME_END is None else min(FRAME_END, total_frames - 1)
if NUM_FRAMES_TO_EXPORT:
    frames_to_export = sorted(set(int(round(x)) for x in np.linspace(FRAME_START, end, NUM_FRAMES_TO_EXPORT)))
else:
    frames_to_export = list(range(FRAME_START, end + 1, FRAME_STRIDE))
print(f"\nExporting {len(frames_to_export)} frames: {frames_to_export[:10]}{' ...' if len(frames_to_export) > 10 else ''}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

rng = np.random.default_rng(RNG_SEED)
frames_meta, value_samples = [], []
center = None
bbox_lo, bbox_hi = None, None

for i, frame_idx in enumerate(frames_to_export):
    data = pipeline.compute(frame_idx)

    pos_all = np.asarray(data.particles.positions, dtype=np.float64)
    vals_all = np.asarray(data.particles[PROPERTY_NAME], dtype=np.float64)
    if vals_all.ndim > 1:
        vals_all = np.linalg.norm(vals_all, axis=1)
    rad_all = get_radii(data, len(pos_all))
    n_raw = len(pos_all)

    if SURFACE_ONLY:
        mask = np.asarray(data.particles['Coordination']) < BULK_COORDINATION
    else:
        mask = np.ones(n_raw, dtype=bool)

    pos, vals, rad = pos_all[mask], vals_all[mask], rad_all[mask]
    ao = compute_ao(pos_all, pos)
    n_vis = len(pos)

    if n_vis > MAX_PARTICLES_PER_FRAME:
        idx = np.sort(rng.choice(n_vis, size=MAX_PARTICLES_PER_FRAME, replace=False))
        pos, vals, rad, ao = pos[idx], vals[idx], rad[idx], ao[idx]
        print(f"  WARNING: frame {frame_idx} has {n_vis:,} surface atoms, subsampled to {MAX_PARTICLES_PER_FRAME:,}")

    if center is None:
        center = pos.mean(axis=0)

    lo, hi = pos.min(axis=0), pos.max(axis=0)
    bbox_lo = lo if bbox_lo is None else np.minimum(bbox_lo, lo)
    bbox_hi = hi if bbox_hi is None else np.maximum(bbox_hi, hi)

    # x, y, z (centred), raw value, radius, ambient occlusion
    out = np.empty((len(pos), 6), dtype=np.float32)
    out[:, 0:3] = pos - center
    out[:, 3] = vals
    out[:, 4] = rad
    out[:, 5] = ao

    bin_name = f"{BASENAME}_frame_{i}.bin"
    out.tofile(os.path.join(OUTPUT_DIR, bin_name))

    frames_meta.append({
        'frame': int(data.attributes.get('Timestep', frame_idx)),
        'file': bin_name, 'count': int(len(pos)),
        'property_min': float(vals.min()), 'property_max': float(vals.max())
    })
    value_samples.append(vals if len(vals) <= 100_000 else rng.choice(vals, 100_000, replace=False))
    print(f"  [{i+1}/{len(frames_to_export)}] frame {frame_idx}: {n_raw:,} atoms -> {n_vis:,} surface -> exported {len(pos):,}")

size = float(np.linalg.norm(bbox_hi - bbox_lo))

samples = np.concatenate(value_samples)
if COLOR_RANGE is not None:
    vmin, vmax = map(float, COLOR_RANGE)
else:
    vmin, vmax = map(float, np.percentile(samples, CLIP_PERCENTILES))
if SYMMETRIC_RANGE:
    m = max(abs(vmin), abs(vmax))
    vmin, vmax = -m, m
print(f"\nColor range ({PROPERTY_NAME}): {vmin:.4g} .. {vmax:.4g}")

meta = {
    'frames': frames_meta, 'size': size,
    'vmin': vmin, 'vmax': vmax,
    'property_name': PROPERTY_NAME,
}
default_cmap = 'jet' if any(k in PROPERTY_NAME.lower() for k in ('pressure', 'stress')) else 'viridis'
bin_total_bytes = sum(os.path.getsize(os.path.join(OUTPUT_DIR, fm['file'])) for fm in frames_meta)

# ----------------------------------------------------------------------------
# HTML template (plain string + placeholders)
# ----------------------------------------------------------------------------
html = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__ - AR Viewer</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: Arial, sans-serif; background: #000; color: #fff; overflow: hidden; }
  #container { width: 100vw; height: 100vh; }
  .panel { position: absolute; background: rgba(0,0,0,0.85); padding: 16px 20px;
           border-radius: 10px; border: 2px solid #0066ff; z-index: 1000; }
  #info { top: 20px; left: 20px; font-size: 14px; line-height: 1.7; max-width: 320px; }
  #info h1 { font-size: 16px; margin-bottom: 8px; }
  #colormap { top: 20px; right: 20px; text-align: center; width: 200px; font-size: 12px; }
  #colormap h3 { font-size: 14px; margin-bottom: 8px; }
  select, input[type=number] { width: 100%; padding: 6px; margin-bottom: 8px;
           background: rgba(255,255,255,0.1); border: 1px solid #0066ff; color: #fff; border-radius: 5px; }
  .cbwrap { display: flex; justify-content: center; gap: 8px; margin-bottom: 10px; }
  #colorbar { width: 40px; height: 180px; border: 2px solid #fff; border-radius: 5px; }
  .cblabels { display: flex; flex-direction: column; justify-content: space-between; font-size: 11px; color: #ccc; text-align: left; }
  .row { text-align: left; color: #aaa; margin-top: 4px; }
  #controls { bottom: 30px; left: 50%; transform: translateX(-50%); text-align: center; }
  button { padding: 10px 18px; margin: 4px; font-size: 14px; background: #0066ff; color: #fff;
           border: none; border-radius: 5px; cursor: pointer; }
  button:hover { background: #0052cc; }
  #frameSlider { width: 400px; margin: 10px auto; }
  #frameInfo { font-size: 16px; font-weight: bold; margin: 8px 0; }
  #loading { position: absolute; top: 50%; left: 50%; transform: translate(-50%,-50%);
             font-size: 24px; background: rgba(0,0,0,0.9); padding: 30px; border-radius: 10px; z-index: 2000; }
  #arNotice { position: absolute; top: 20px; left: 50%; transform: translateX(-50%);
              background: rgba(0,0,0,0.85); padding: 10px 18px; border-radius: 8px;
              font-size: 12px; color: #7fdcff; z-index: 999; text-align: center; }
</style>
</head>
<body>
<div id="loading">Loading particle data...</div>
<div id="container"></div>

<div id="info" class="panel" style="display:none;">
  <h1>__TITLE__</h1>
  <div id="particleInfo">Loading...</div>
  <div style="font-size:11px; color:#aaa; margin-top:10px; line-height:1.5;">
    <b>Hand controls (in AR):</b><br>
    Turn a palm face-up: model rests in your hand<br>
    Pinch inside the model + move one hand: reposition<br>
    Pinch with both hands: pull apart/together to scale, twist to rotate
  </div>
</div>

<div id="colormap" class="panel" style="display:none;">
  <h3>Colormap</h3>
  <select id="cmap" onchange="updateColormap()">
    <option value="viridis">Viridis</option>
    <option value="plasma">Plasma</option>
    <option value="inferno">Inferno</option>
    <option value="magma">Magma</option>
    <option value="jet">Jet</option>
    <option value="hot">Hot</option>
    <option value="cool">Cool</option>
  </select>
  <div class="cbwrap">
    <div id="colorbar"></div>
    <div class="cblabels"><span id="cbMax"></span><span id="cbMid"></span><span id="cbMin"></span></div>
  </div>
  <div class="row">Max</div><input id="rangeMax" type="number" step="any" onchange="updateRange()">
  <div class="row">Min</div><input id="rangeMin" type="number" step="any" onchange="updateRange()">
  <button onclick="resetRange()">Auto range</button>
  <div class="row">Atom size</div>
  <input id="sizeSlider" type="range" min="0.3" max="2" step="0.05" value="1" oninput="updateAtomSize()" style="width:100%">
</div>

<div id="controls" class="panel" style="display:none;">
  <div id="frameInfo">Frame 1 / __NFRAMES__</div>
  <div>
    <button onclick="firstFrame()">First</button>
    <button onclick="prevFrame()">Previous</button>
    <button id="playBtn" onclick="togglePlay()">Play</button>
    <button onclick="nextFrame()">Next</button>
    <button onclick="lastFrame()">Last</button>
    <button onclick="recenterModel()">Recenter</button>
    <button id="palmModeBtn" onclick="togglePalmMode()">Palm Mode: ON</button>
    <button id="upAxisBtn" onclick="cycleUpAxis()">Up axis</button>
  </div>
  <input type="range" id="frameSlider" min="0" max="__NFRAMES_MINUS_1__" value="0" oninput="onSliderChange()">
  <div style="font-size:12px; margin-top:8px;">Desktop preview: drag to rotate, scroll to zoom, arrow keys step, space plays</div>
</div>

<div id="arNotice" style="display:none;">
  Enable Hand Tracking in Quest Settings &gt; Movement Tracking, then tap "Enter AR" below
</div>

<script type="importmap">
{"imports": {"three": "https://unpkg.com/three@0.160.0/build/three.module.js",
             "three/addons/": "https://unpkg.com/three@0.160.0/examples/jsm/"}}
</script>

<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { ARButton } from 'three/addons/webxr/ARButton.js';

const META = __META_JSON__;
const STRIDE = 6;              // x, y, z, value, radius, ao
const SHOW_HELPERS = false;    // grid + axes
const PLAY_INTERVAL_MS = __PLAY_INTERVAL_MS__;
const MAX_CACHED_FRAMES = __MAX_CACHED_FRAMES__;
const PREFETCH_AHEAD = __PREFETCH_AHEAD__;

// ---------------- colormaps ----------------
const COLORMAPS = {
  viridis: [[0.267,0.005,0.329],[0.283,0.141,0.458],[0.254,0.265,0.530],[0.164,0.471,0.558],[0.135,0.659,0.518],[0.478,0.821,0.318],[0.993,0.906,0.144]],
  plasma:  [[0.050,0.030,0.528],[0.417,0.001,0.658],[0.693,0.165,0.565],[0.881,0.393,0.383],[0.987,0.810,0.145]],
  inferno: [[0.001,0.000,0.014],[0.144,0.047,0.338],[0.502,0.092,0.472],[0.855,0.258,0.283],[0.988,0.553,0.111]],
  magma:   [[0.001,0.000,0.014],[0.129,0.046,0.285],[0.429,0.101,0.478],[0.670,0.188,0.427],[0.892,0.364,0.279],[0.987,0.810,0.145]],
  jet:  t => { const r=Math.max(0,Math.min(1,1.5-Math.abs(4*t-3))), g=Math.max(0,Math.min(1,1.5-Math.abs(4*t-2))), b=Math.max(0,Math.min(1,1.5-Math.abs(4*t-1))); return [r,g,b]; },
  hot:  t => [Math.min(1,3*t), Math.min(1,Math.max(0,3*t-1)), Math.min(1,Math.max(0,3*t-2))],
  cool: t => [t, 1-t, 1]
};
function getColor(cmap, t) {
  if (typeof cmap === 'function') return cmap(t);
  const i = Math.min(Math.floor(t*(cmap.length-1)), cmap.length-2), f = t*(cmap.length-1) - i;
  return cmap[i].map((v,k) => v + (cmap[i+1][k] - v) * f);
}

// 256-entry lookup texture -> colormap changes are instant (no reload)
const lutData = new Uint8Array(256 * 4);
const lutTex = new THREE.DataTexture(lutData, 256, 1, THREE.RGBAFormat);
lutTex.magFilter = THREE.LinearFilter;
lutTex.minFilter = THREE.LinearFilter;
lutTex.wrapS = lutTex.wrapT = THREE.ClampToEdgeWrapping;
function fillLUT(name) {
  const cm = COLORMAPS[name];
  for (let i = 0; i < 256; i++) {
    const c = getColor(cm, i / 255);
    lutData[i*4] = Math.round(c[0]*255); lutData[i*4+1] = Math.round(c[1]*255);
    lutData[i*4+2] = Math.round(c[2]*255); lutData[i*4+3] = 255;
  }
  lutTex.needsUpdate = true;
}

// ---------------- sphere impostor shader ----------------
const vertexShader = `
  attribute float aValue;
  attribute float aRadius;
  attribute float aAO;
  uniform float uViewportHeight;
  uniform float uRadiusScale;
  varying float vValue;
  varying float vAO;
  varying vec3  vViewCenter;
  varying float vViewRadius;
  void main() {
    vec4 mv = modelViewMatrix * vec4(position, 1.0);
    float worldScale = length(modelMatrix[0].xyz);
    float r = aRadius * uRadiusScale * worldScale;
    vViewCenter = mv.xyz;
    vViewRadius = r;
    vValue = aValue;
    vAO = aAO;
    gl_Position = projectionMatrix * mv;
    gl_PointSize = max(1.0, uViewportHeight * projectionMatrix[1][1] * r / max(-mv.z, 1e-4));
  }`;

const fragmentShader = `
  uniform sampler2D uLUT;
  uniform float uMin;
  uniform float uMax;
  uniform vec3  uLightDir;
  uniform mat4  projectionMatrix;
  varying float vValue;
  varying float vAO;
  varying vec3  vViewCenter;
  varying float vViewRadius;
  void main() {
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    c.y = -c.y;
    float r2 = dot(c, c);
    if (r2 > 1.0) discard;
    vec3 n = vec3(c, sqrt(1.0 - r2));

    // correct per-pixel depth -> spheres intersect like real spheres
    vec3 p = vViewCenter + n * vViewRadius;
    vec4 clip = projectionMatrix * vec4(p, 1.0);
    gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;

    float t = clamp((vValue - uMin) / max(uMax - uMin, 1e-20), 0.0, 1.0);
    vec3 base = texture2D(uLUT, vec2(t, 0.5)).rgb;

    float diff = max(dot(n, uLightDir), 0.0);
    vec3  h    = normalize(uLightDir + vec3(0.0, 0.0, 1.0));
    float spec = pow(max(dot(n, h), 0.0), 48.0);
    float rim  = mix(0.55, 1.0, n.z);     // darker silhouette, like Tachyon

    vec3 col = base * (0.30 + 0.75 * diff) * vAO * rim + vec3(0.35) * spec * vAO;
    gl_FragColor = vec4(col, 1.0);
  }`;

const sphereMat = new THREE.ShaderMaterial({
  uniforms: {
    uLUT:            { value: lutTex },
    uMin:            { value: META.vmin },
    uMax:            { value: META.vmax },
    uViewportHeight: { value: 1 },
    uRadiusScale:    { value: 1 },
    uLightDir:       { value: new THREE.Vector3(-0.35, 0.55, 1.0).normalize() }
  },
  vertexShader, fragmentShader
});

// ---------------- scene ----------------
let particleSystem = null;
let currentCmap = '__DEFAULT_CMAP__', currentFrameIdx = 0, isPlaying = false;

const DESIRED_WORLD_SIZE = 0.35;
const DEFAULT_POSITION = new THREE.Vector3(0, 1.2, -0.45);
const initialScale = DESIRED_WORLD_SIZE / META.size;
const baseHalfExtent = META.size * 0.5;

const scene = new THREE.Scene();
scene.background = null;
const camera = new THREE.PerspectiveCamera(60, window.innerWidth/window.innerHeight, 0.01, 100);

const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0x000000, 0);
renderer.xr.enabled = true;
renderer.xr.setReferenceSpaceType('local-floor');
renderer.xr.setFoveation(0);   // full resolution everywhere (sharper atoms)
document.getElementById('container').appendChild(renderer.domElement);

document.body.appendChild(ARButton.createButton(renderer, {
  optionalFeatures: ['hand-tracking', 'local-floor', 'dom-overlay'],
  domOverlay: { root: document.body }
}));

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.05;
camera.position.set(0, 1.2, 0.6);
controls.target.copy(DEFAULT_POSITION);
controls.update();

const modelGroup = new THREE.Group();
modelGroup.position.copy(DEFAULT_POSITION);
modelGroup.scale.setScalar(initialScale);
scene.add(modelGroup);

// Inner group that turns the simulation's 'up' axis into the AR world's up (+Y).
// Hand rotation/recenter act on modelGroup, so this orientation is never lost.
const orientGroup = new THREE.Group();
modelGroup.add(orientGroup);
const UP_AXES = ['+z', '-z', '+x', '-x', '+y', '-y'];
let upAxisIdx = Math.max(0, UP_AXES.indexOf('__UP_AXIS__'));
function applyUpAxis() {
  const a = UP_AXES[upAxisIdx];
  orientGroup.rotation.set(0, 0, 0);
  if (a === '+z') orientGroup.rotation.x = -Math.PI / 2;
  else if (a === '-z') orientGroup.rotation.x = Math.PI / 2;
  else if (a === '+x') orientGroup.rotation.z = Math.PI / 2;
  else if (a === '-x') orientGroup.rotation.z = -Math.PI / 2;
  else if (a === '-y') orientGroup.rotation.x = Math.PI;
  document.getElementById('upAxisBtn').textContent = `Up axis: ${a.toUpperCase()}`;
}
window.cycleUpAxis = () => { upAxisIdx = (upAxisIdx + 1) % UP_AXES.length; applyUpAxis(); };

if (SHOW_HELPERS) {
  const refGrid = new THREE.GridHelper(META.size * 1.2, 10, 0x0066ff, 0x003366);
  refGrid.position.y = -META.size * 0.5;
  modelGroup.add(refGrid);
  modelGroup.add(new THREE.AxesHelper(META.size * 0.3));
}

// ---------------- frame loading (LRU cache + prefetch) ----------------
const frameCache = new Map();

function evictCache() {
  let guard = frameCache.size;
  while (frameCache.size > MAX_CACHED_FRAMES && guard-- > 0) {
    const oldest = frameCache.keys().next().value;
    const p = frameCache.get(oldest);
    frameCache.delete(oldest);
    if (oldest === currentFrameIdx) { frameCache.set(oldest, p); continue; }
    p.then(pts => { if (pts !== particleSystem) pts.geometry.dispose(); }).catch(() => {});
  }
}

function loadFrame(i) {
  if (frameCache.has(i)) {             // mark as recently used
    const p = frameCache.get(i);
    frameCache.delete(i); frameCache.set(i, p);
    return p;
  }
  const p = (async () => {
    const m = META.frames[i];
    const resp = await fetch(m.file);
    if (!resp.ok) throw new Error(`Could not load ${m.file} (HTTP ${resp.status})`);
    const buf = await resp.arrayBuffer();
    const ib = new THREE.InterleavedBuffer(new Float32Array(buf), STRIDE);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.InterleavedBufferAttribute(ib, 3, 0));
    geo.setAttribute('aValue',   new THREE.InterleavedBufferAttribute(ib, 1, 3));
    geo.setAttribute('aRadius',  new THREE.InterleavedBufferAttribute(ib, 1, 4));
    geo.setAttribute('aAO',      new THREE.InterleavedBufferAttribute(ib, 1, 5));
    const pts = new THREE.Points(geo, sphereMat);
    pts.frustumCulled = false;
    return pts;
  })();
  frameCache.set(i, p);
  p.catch(() => frameCache.delete(i));
  evictCache();
  return p;
}

function prefetch(idx) {
  for (let k = 1; k <= PREFETCH_AHEAD; k++) loadFrame((idx + k) % META.frames.length);
}

function fmt(v) {
  const a = Math.abs(v);
  return (a !== 0 && (a >= 1e4 || a < 1e-2)) ? v.toExponential(2) : v.toFixed(2);
}

let frameLoading = false;
async function updateFrame() {
  const idx = currentFrameIdx;
  const meta = META.frames[idx];
  if (!frameCache.has(idx))
    document.getElementById('particleInfo').innerHTML = `Loading ${meta.count.toLocaleString()} atoms...`;
  frameLoading = true;
  let pts;
  try { pts = await loadFrame(idx); }
  catch (err) {
    frameLoading = false;
    document.getElementById('particleInfo').innerHTML = `<span style="color:#ff6b6b">${err.message}</span>`;
    return;
  }
  frameLoading = false;
  if (idx !== currentFrameIdx) return;   // a newer frame was requested meanwhile
  if (particleSystem) orientGroup.remove(particleSystem);
  particleSystem = pts;
  orientGroup.add(pts);
  document.getElementById('particleInfo').innerHTML =
    `Atoms shown: ${meta.count.toLocaleString()}<br>Property: ${META.property_name}<br>` +
    `Data range: ${fmt(meta.property_min)} .. ${fmt(meta.property_max)}`;
  const label = `Frame ${idx+1} / ${META.frames.length} (Step ${meta.frame})`;
  document.getElementById('frameInfo').textContent = label;
  document.getElementById('frameSlider').value = idx;
  if (typeof arFrameLabel !== 'undefined') arFrameLabel.setText(label);
  prefetch(idx);
}

function updateGradient() {
  const cmap = COLORMAPS[currentCmap];
  let grad = 'linear-gradient(to bottom';
  for (let i = 0; i <= 100; i++) {
    const c = getColor(cmap, 1 - i/100);
    grad += `,rgb(${Math.round(c[0]*255)},${Math.round(c[1]*255)},${Math.round(c[2]*255)}) ${i}%`;
  }
  document.getElementById('colorbar').style.background = grad + ')';
  const u = sphereMat.uniforms;
  document.getElementById('cbMax').textContent = fmt(u.uMax.value);
  document.getElementById('cbMid').textContent = fmt(0.5*(u.uMin.value + u.uMax.value));
  document.getElementById('cbMin').textContent = fmt(u.uMin.value);
  document.getElementById('rangeMax').value = u.uMax.value;
  document.getElementById('rangeMin').value = u.uMin.value;
}

window.updateColormap = () => { currentCmap = document.getElementById('cmap').value; fillLUT(currentCmap); updateGradient(); };
window.updateRange = () => {
  const lo = parseFloat(document.getElementById('rangeMin').value);
  const hi = parseFloat(document.getElementById('rangeMax').value);
  if (isFinite(lo) && isFinite(hi) && hi > lo) { sphereMat.uniforms.uMin.value = lo; sphereMat.uniforms.uMax.value = hi; }
  updateGradient();
};
window.resetRange = () => { sphereMat.uniforms.uMin.value = META.vmin; sphereMat.uniforms.uMax.value = META.vmax; updateGradient(); };
window.updateAtomSize = () => { sphereMat.uniforms.uRadiusScale.value = parseFloat(document.getElementById('sizeSlider').value); };

window.firstFrame = () => { currentFrameIdx = 0; updateFrame(); };
window.lastFrame  = () => { currentFrameIdx = META.frames.length-1; updateFrame(); };
window.nextFrame  = () => { currentFrameIdx = (currentFrameIdx+1) % META.frames.length; updateFrame(); };
window.prevFrame  = () => { currentFrameIdx = (currentFrameIdx-1+META.frames.length) % META.frames.length; updateFrame(); };
window.onSliderChange = () => { currentFrameIdx = parseInt(document.getElementById('frameSlider').value); updateFrame(); };
window.recenterModel = () => {
  modelGroup.position.copy(DEFAULT_POSITION);
  modelGroup.scale.setScalar(initialScale);
  modelGroup.rotation.set(0,0,0);
};
let lastAdvance = 0;
window.togglePlay = () => {
  isPlaying = !isPlaying;
  document.getElementById('playBtn').textContent = isPlaying ? 'Pause' : 'Play';
  if (typeof arPlayBtn !== 'undefined') arPlayBtn.setText(isPlaying ? 'Pause' : 'Play');
};
document.addEventListener('keydown', e => {
  if (e.key === 'ArrowLeft') window.prevFrame();
  if (e.key === 'ArrowRight') window.nextFrame();
  if (e.key === ' ') { e.preventDefault(); window.togglePlay(); }
});

// ---------------- hand tracking ----------------
const hand1 = renderer.xr.getHand(0);
const hand2 = renderer.xr.getHand(1);
scene.add(hand1); scene.add(hand2);

function makePinchMarker() {
  const m = new THREE.Mesh(new THREE.SphereGeometry(0.008, 12, 12),
    new THREE.MeshBasicMaterial({ color: 0x00ffff, transparent: true, opacity: 0.85 }));
  m.visible = false; scene.add(m); return m;
}
const marker1 = makePinchMarker();
const marker2 = makePinchMarker();

const PINCH_THRESHOLD = 0.028;
const _tmpThumb = new THREE.Vector3(), _tmpIndex = new THREE.Vector3();
function getPinchPoint(hand, out) {
  const thumb = hand.joints['thumb-tip'], index = hand.joints['index-finger-tip'];
  if (!thumb || !index) return null;
  thumb.getWorldPosition(_tmpThumb); index.getWorldPosition(_tmpIndex);
  if (_tmpThumb.distanceTo(_tmpIndex) > PINCH_THRESHOLD) return null;
  return out.addVectors(_tmpThumb, _tmpIndex).multiplyScalar(0.5);
}

const pinch1 = new THREE.Vector3(), pinch2 = new THREE.Vector3();
let singleGrab = null, twoHandGrab = null;
const worldPos = new THREE.Vector3();

let hand1Handedness = null, hand2Handedness = null;
hand1.addEventListener('connected', e => { hand1Handedness = e.data.handedness; });
hand1.addEventListener('disconnected', () => { hand1Handedness = null; });
hand2.addEventListener('connected', e => { hand2Handedness = e.data.handedness; });
hand2.addEventListener('disconnected', () => { hand2Handedness = null; });

let palmModeEnabled = true;
window.togglePalmMode = () => {
  palmModeEnabled = !palmModeEnabled;
  document.getElementById('palmModeBtn').textContent = `Palm Mode: ${palmModeEnabled ? 'ON' : 'OFF'}`;
};

const WORLD_UP = new THREE.Vector3(0,1,0);
const PALM_UP_THRESHOLD = 0.55, PALM_HOVER_OFFSET = 0.05, PALM_FOLLOW_LERP = 0.18;
const PALM_NORMAL_SIGN = { left: 1, right: -1 };
const _palmWrist = new THREE.Vector3(), _palmIndex = new THREE.Vector3(), _palmPinky = new THREE.Vector3();
const _palmPos = new THREE.Vector3(), _palmNormal = new THREE.Vector3();
const _palmDirA = new THREE.Vector3(), _palmDirB = new THREE.Vector3();

function getPalmPose(hand, handedness) {
  const wrist = hand.joints['wrist'], indexMeta = hand.joints['index-finger-metacarpal'], pinkyMeta = hand.joints['pinky-finger-metacarpal'];
  if (!wrist || !indexMeta || !pinkyMeta || !handedness) return null;
  wrist.getWorldPosition(_palmWrist); indexMeta.getWorldPosition(_palmIndex); pinkyMeta.getWorldPosition(_palmPinky);
  _palmPos.addVectors(_palmIndex, _palmPinky).multiplyScalar(0.5).add(_palmWrist).multiplyScalar(0.5);
  _palmDirA.subVectors(_palmIndex, _palmWrist); _palmDirB.subVectors(_palmPinky, _palmWrist);
  _palmNormal.crossVectors(_palmDirA, _palmDirB).normalize().multiplyScalar(PALM_NORMAL_SIGN[handedness] || 1);
  return { position: _palmPos, normal: _palmNormal };
}

function updatePalmFollow() {
  if (!palmModeEnabled || singleGrab || twoHandGrab) return false;
  for (const [hand, handedness] of [[hand1, hand1Handedness], [hand2, hand2Handedness]]) {
    const pose = getPalmPose(hand, handedness);
    if (pose && pose.normal.dot(WORLD_UP) > PALM_UP_THRESHOLD) {
      const target = pose.position.clone().addScaledVector(pose.normal, PALM_HOVER_OFFSET);
      modelGroup.position.lerp(target, PALM_FOLLOW_LERP);
      return true;
    }
  }
  return false;
}

function updateHandInteraction() {
  if (!renderer.xr.isPresenting) return;
  const onPalm = updatePalmFollow();
  const p1 = getPinchPoint(hand1, pinch1), p2 = getPinchPoint(hand2, pinch2);
  marker1.visible = !!p1; if (p1) marker1.position.copy(p1);
  marker2.visible = !!p2; if (p2) marker2.position.copy(p2);
  modelGroup.getWorldPosition(worldPos);
  const worldHalfExtent = baseHalfExtent * modelGroup.scale.x;
  if (onPalm) return;

  if (p1 && p2) {
    const dist = p1.distanceTo(p2);
    const mid = new THREE.Vector3().addVectors(p1, p2).multiplyScalar(0.5);
    const angle = Math.atan2(p2.z - p1.z, p2.x - p1.x);
    const nearModel = Math.min(p1.distanceTo(worldPos), p2.distanceTo(worldPos)) < worldHalfExtent * 1.6;
    if (!twoHandGrab && nearModel) {
      twoHandGrab = { startDist: dist, startAngle: angle, startScale: modelGroup.scale.x,
                      startRotY: modelGroup.rotation.y, startMid: mid.clone(), startPos: modelGroup.position.clone() };
    } else if (twoHandGrab) {
      const s = THREE.MathUtils.clamp(twoHandGrab.startScale * (dist / twoHandGrab.startDist), initialScale*0.15, initialScale*6);
      modelGroup.scale.setScalar(s);
      modelGroup.rotation.y = twoHandGrab.startRotY + (angle - twoHandGrab.startAngle);
      modelGroup.position.copy(twoHandGrab.startPos).add(new THREE.Vector3().subVectors(mid, twoHandGrab.startMid));
    }
    singleGrab = null;
  } else {
    twoHandGrab = null;
    const active = p1 ? { which: 1, pt: p1 } : (p2 ? { which: 2, pt: p2 } : null);
    if (active) {
      if (!singleGrab && active.pt.distanceTo(worldPos) < worldHalfExtent * 1.4) {
        singleGrab = { which: active.which, offset: new THREE.Vector3().subVectors(modelGroup.position, active.pt) };
      } else if (singleGrab && singleGrab.which === active.which) {
        modelGroup.position.copy(active.pt).add(singleGrab.offset);
      }
    } else {
      singleGrab = null;
    }
  }
}

// ---------------- in-AR button panel (poke with your index finger) ----------------
const arPanel = new THREE.Group();
arPanel.visible = false;
scene.add(arPanel);
const arButtons = [];

function drawLabel(canvas, text, bg) {
  const g = canvas.getContext('2d');
  g.clearRect(0, 0, canvas.width, canvas.height);
  g.fillStyle = bg;
  g.fillRect(0, 0, canvas.width, canvas.height);
  g.strokeStyle = '#ffffff';
  g.lineWidth = 6;
  g.strokeRect(3, 3, canvas.width - 6, canvas.height - 6);
  g.fillStyle = '#ffffff';
  g.font = `bold ${Math.round(canvas.height * 0.42)}px Arial`;
  g.textAlign = 'center';
  g.textBaseline = 'middle';
  g.fillText(text, canvas.width / 2, canvas.height / 2);
}

function makePanelItem(text, w, h, x, y, bg, onPress) {
  const canvas = document.createElement('canvas');
  canvas.width = Math.round(w * 4000);
  canvas.height = Math.round(h * 4000);
  drawLabel(canvas, text, bg);
  const tex = new THREE.CanvasTexture(canvas);
  tex.colorSpace = THREE.SRGBColorSpace;
  const mesh = new THREE.Mesh(new THREE.PlaneGeometry(w, h),
    new THREE.MeshBasicMaterial({ map: tex, transparent: true, opacity: 0.92, depthTest: false }));
  mesh.position.set(x, y, 0);
  mesh.renderOrder = 999;
  arPanel.add(mesh);
  const item = { mesh, w, h, onPress, pressed: false,
                 setText(t) { drawLabel(canvas, t, bg); tex.needsUpdate = true; } };
  if (onPress) arButtons.push(item);
  return item;
}

const BW = 0.052, BH = 0.036, GAP = 0.058;
const arFrameLabel = makePanelItem('Frame 1', 0.29, 0.03, 0, 0.042, '#000000', null);
makePanelItem('<<',    BW, BH, -2*GAP, 0, '#0066ff', () => window.prevFrame());
const arPlayBtn =
makePanelItem('Play',  BW, BH, -1*GAP, 0, '#0066ff', () => window.togglePlay());
makePanelItem('>>',    BW, BH,  0,     0, '#0066ff', () => window.nextFrame());
makePanelItem('Up',    BW, BH,  1*GAP, 0, '#444444', () => window.cycleUpAxis());
makePanelItem('Reset', BW, BH,  2*GAP, 0, '#444444', () => window.recenterModel());

const _head = new THREE.Vector3(), _fwd = new THREE.Vector3();
const _tip = new THREE.Vector3(), _local = new THREE.Vector3();
let panelPlaced = false, lastPokeTime = 0, xrFrames = 0;

function placePanel() {
  camera.getWorldPosition(_head);
  camera.getWorldDirection(_fwd);
  _fwd.y = 0;
  if (_fwd.lengthSq() < 1e-6) _fwd.set(0, 0, -1);
  _fwd.normalize();
  arPanel.position.copy(_head).addScaledVector(_fwd, 0.38);
  arPanel.position.y -= 0.28;                       // chest height, below your view
  arPanel.lookAt(_head.x, arPanel.position.y, _head.z);
  panelPlaced = true;
}

function updateArPanel(time) {
  arPanel.visible = renderer.xr.isPresenting;
  if (!arPanel.visible) { panelPlaced = false; xrFrames = 0; return; }
  if (!panelPlaced) { if (++xrFrames > 5) placePanel(); return; }   // wait until the head pose is known
  camera.getWorldPosition(_head);
  if (_head.distanceTo(arPanel.position) > 1.5) placePanel();   // follow if you walk away

  for (const b of arButtons) {
    let touching = false;
    for (const hand of [hand1, hand2]) {
      const tip = hand.joints['index-finger-tip'];
      if (!tip || !tip.visible) continue;
      tip.getWorldPosition(_tip);
      _local.copy(_tip);
      b.mesh.worldToLocal(_local);
      if (Math.abs(_local.x) < b.w / 2 && Math.abs(_local.y) < b.h / 2 && Math.abs(_local.z) < 0.02) { touching = true; break; }
    }
    if (touching && !b.pressed && time - lastPokeTime > 400) {
      b.pressed = true; lastPokeTime = time;
      b.mesh.scale.setScalar(0.85);
      b.onPress();
    } else if (!touching && b.pressed) {
      b.pressed = false;
      b.mesh.scale.setScalar(1);
    }
  }
}

renderer.xr.addEventListener('sessionstart', () => { controls.enabled = false; });
renderer.xr.addEventListener('sessionend', () => { controls.enabled = true; singleGrab = null; twoHandGrab = null; });

// point size must know the pixel height of the current view (desktop or per-eye in XR)
function updateViewportHeight() {
  let h = renderer.domElement.height;
  if (renderer.xr.isPresenting) {
    const xc = renderer.xr.getCamera();
    if (xc.cameras.length && xc.cameras[0].viewport.w > 0) h = xc.cameras[0].viewport.w;
  }
  sphereMat.uniforms.uViewportHeight.value = h;
}

async function init() {
  document.getElementById('cmap').value = currentCmap;
  fillLUT(currentCmap);
  applyUpAxis();
  await updateFrame();
  updateGradient();
  for (const id of ['info', 'colormap', 'controls', 'arNotice'])
    document.getElementById(id).style.display = 'block';
  document.getElementById('loading').style.display = 'none';
}
init();

renderer.setAnimationLoop((time) => {
  if (isPlaying && !frameLoading && time - lastAdvance >= PLAY_INTERVAL_MS) {
    lastAdvance = time;
    window.nextFrame();
  }
  updateHandInteraction();
  updateArPanel(time);
  if (controls.enabled) controls.update();
  updateViewportHeight();
  renderer.render(scene, camera);
});
window.addEventListener('resize', () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});
</script>
</body>
</html>"""

html = (html.replace("__META_JSON__", json.dumps(meta))
            .replace("__DEFAULT_CMAP__", default_cmap)
            .replace("__UP_AXIS__", UP_AXIS)
            .replace("__PLAY_INTERVAL_MS__", str(int(PLAY_INTERVAL_MS)))
            .replace("__MAX_CACHED_FRAMES__", str(int(MAX_CACHED_FRAMES)))
            .replace("__PREFETCH_AHEAD__", str(int(PREFETCH_AHEAD)))
            .replace("__NFRAMES_MINUS_1__", str(len(frames_meta) - 1))
            .replace("__NFRAMES__", str(len(frames_meta)))
            .replace("__TITLE__", TITLE))

html_file = os.path.join(OUTPUT_DIR, f'{BASENAME}_ar.html')
with open(html_file, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n{'='*70}")
print("AR VIEWER CREATED")
print(f"{'='*70}")
print(f"Saved: {html_file}")
print(f"Plus {len(frames_meta)} .bin files ({bin_total_bytes/1024/1024:.1f} MB total)")
print("\nPreview on PC:")
print(f'  cd "{OUTPUT_DIR}"')
print("  python -m http.server 8000")
print(f"  open http://localhost:8000/{BASENAME}_ar.html")
print(f"{'='*70}")
