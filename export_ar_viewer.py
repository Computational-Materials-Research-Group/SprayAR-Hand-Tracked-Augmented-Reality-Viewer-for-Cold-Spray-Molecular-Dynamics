# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 17:16:12 2026

@author: pedit
"""

import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI.*')
import ovito._extensions.pyscript

from ovito.io import import_file
import numpy as np
import os
import json

print("="*70)
print("COLD SPRAY DEPOSITION - AR HAND-TRACKING VIEWER EXPORT")
print("="*70)

INPUT_FILE = r"D:\Atomisitc modeling research\Electrostatic cold spray deposition\deposition_postimpact.lammpstrj"
OUTPUT_DIR = os.path.join(os.path.dirname(INPUT_FILE), "ar_output")
BASENAME = "deposition_postimpact"

PROPERTY_NAME = "v_pressure"

MAX_PARTICLES_PER_FRAME = 200_000
RNG_SEED = 42

NUM_FRAMES_TO_EXPORT = 6

pipeline = import_file(INPUT_FILE)
total_frames = pipeline.source.num_frames
print(f"Total frames: {total_frames}")

data0 = pipeline.compute(0)
available_props = list(data0.particles.keys())
print(f"Available particle properties: {available_props}")

if PROPERTY_NAME not in available_props:
    fallback = "Particle Type" if "Particle Type" in available_props else available_props[0]
    print(f"Property '{PROPERTY_NAME}' not found in this file -> using '{fallback}' instead.")
    print("(Set PROPERTY_NAME to one of the properties above and re-run to change this.)")
    PROPERTY_NAME = fallback

frames_to_export = sorted(set(int(round(x)) for x in np.linspace(0, total_frames - 1, NUM_FRAMES_TO_EXPORT)))
print(f"\nExporting {len(frames_to_export)} frames (frame indices: {frames_to_export}), "
      f"capped at {MAX_PARTICLES_PER_FRAME:,} particles/frame...")

output_dir = OUTPUT_DIR
os.makedirs(output_dir, exist_ok=True)

frames_meta = []
global_min, global_max = None, None
center, size = None, None

for i, frame_idx in enumerate(frames_to_export):
    data = pipeline.compute(frame_idx)

    positions = np.array(data.particles.positions)
    property_values = np.array(data.particles[PROPERTY_NAME])
    if property_values.ndim > 1:
        property_values = np.linalg.norm(property_values, axis=1)
    property_values = property_values.astype(np.float64)

    n_raw = len(positions)
    if n_raw > MAX_PARTICLES_PER_FRAME:
        rng = np.random.default_rng(RNG_SEED)
        idx = np.sort(rng.choice(n_raw, size=MAX_PARTICLES_PER_FRAME, replace=False))
        positions = positions[idx]
        property_values = property_values[idx]

    pmin, pmax = float(property_values.min()), float(property_values.max())
    prop_norm = (property_values - pmin) / (pmax - pmin) if pmax > pmin else np.ones_like(property_values) * 0.5

    if i == 0:
        center = positions.mean(axis=0)
        size = float(np.linalg.norm(positions.max(axis=0) - positions.min(axis=0)))

    interleaved = np.empty((len(positions), 4), dtype=np.float32)
    interleaved[:, 0:3] = positions
    interleaved[:, 3] = prop_norm

    bin_name = f"{BASENAME}_frame_{i}.bin"
    interleaved.tofile(os.path.join(output_dir, bin_name))

    frames_meta.append({
        'frame': frame_idx,
        'file': bin_name,
        'count': int(len(positions)),
        'property_min': pmin,
        'property_max': pmax
    })

    global_min = pmin if global_min is None else min(global_min, pmin)
    global_max = pmax if global_max is None else max(global_max, pmax)

    print(f"  Frame {frame_idx}: {n_raw:,} particles in trajectory -> exported {len(positions):,} -> {bin_name}")

meta = {
    'frames': frames_meta,
    'center': center.tolist(),
    'size': size,
    'global_min': global_min,
    'global_max': global_max,
    'property_name': PROPERTY_NAME,
}
meta_json = json.dumps(meta)

default_cmap = 'jet' if any(k in PROPERTY_NAME.lower() for k in ('pressure', 'stress')) else 'viridis'
print(f"Default colormap: {default_cmap} (change anytime from the in-viewer dropdown)")

bin_total_bytes = sum(os.path.getsize(os.path.join(output_dir, fm['file'])) for fm in frames_meta)
print(f"\nTotal binary data size: {bin_total_bytes/1024/1024:.1f} MB across {len(frames_meta)} files")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Cold Spray Deposition - Post-Impact - AR Viewer</title>
    <style>
        * {{ margin: 0; padding: 0; }}
        body {{ font-family: Arial, sans-serif; background: #000; color: #fff; overflow: hidden; }}
        #container {{ width: 100vw; height: 100vh; }}
        #info {{
            position: absolute; top: 20px; left: 20px;
            background: rgba(0,0,0,0.85); padding: 20px; border-radius: 10px;
            border: 2px solid #0066ff; z-index: 1000; font-size: 14px; line-height: 1.8;
            max-width: 320px;
        }}
        #info h1 {{ font-size: 16px; margin-bottom: 10px; }}
        #colormap {{
            position: absolute; top: 20px; right: 20px;
            background: rgba(0,0,0,0.85); padding: 20px; border-radius: 10px;
            border: 2px solid #0066ff; z-index: 1000; text-align: center;
        }}
        select {{
            width: 150px; padding: 8px; margin-bottom: 15px;
            background: rgba(255,255,255,0.1); border: 1px solid #0066ff;
            color: white; border-radius: 5px; cursor: pointer;
        }}
        #colorbar {{ width: 50px; height: 200px; border: 2px solid #fff; border-radius: 5px; margin: 0 auto 10px; }}
        .label {{ font-size: 11px; color: #aaa; }}
        #controls {{
            position: absolute; bottom: 30px; left: 50%; transform: translateX(-50%);
            background: rgba(0,0,0,0.85); padding: 20px 30px; border-radius: 10px;
            border: 2px solid #0066ff; z-index: 1000; text-align: center;
        }}
        button {{
            padding: 10px 20px; margin: 5px; font-size: 14px;
            background: #0066ff; color: white; border: none;
            border-radius: 5px; cursor: pointer;
        }}
        button:hover {{ background: #0052cc; }}
        #frameSlider {{ width: 400px; margin: 10px auto; }}
        #frameInfo {{ font-size: 16px; font-weight: bold; margin: 10px 0; }}
        #loading {{
            position: absolute; top: 50%; left: 50%; transform: translate(-50%,-50%);
            font-size: 24px; background: rgba(0,0,0,0.9); padding: 30px;
            border-radius: 10px; z-index: 2000;
        }}
        #arNotice {{
            position: absolute; bottom: 140px; left: 50%; transform: translateX(-50%);
            background: rgba(0,0,0,0.85); padding: 10px 18px; border-radius: 8px;
            font-size: 12px; color: #7fdcff; z-index: 999; text-align: center;
        }}
    </style>
</head>
<body>
    <div id="loading">Loading particle data...</div>
    <div id="container"></div>

    <div id="info" style="display:none;">
        <h1>Cold Spray Deposition - Post-Impact</h1>
        <div id="particleInfo">Loading...</div>
        <div style="font-size:11px; color:#aaa; margin-top:10px; line-height:1.5;">
            <b>Hand controls (in AR):</b><br>
            Turn a palm face-up: model rests in your hand<br>
            Pinch inside the model + move one hand: reposition<br>
            Pinch with both hands: pull apart/together to scale, twist to rotate
        </div>
    </div>

    <div id="colormap" style="display:none;">
        <h3 style="font-size:14px; margin-bottom:10px;">Colormap</h3>
        <select id="cmap" onchange="updateColormap()">
            <option value="viridis">Viridis</option>
            <option value="plasma">Plasma</option>
            <option value="inferno">Inferno</option>
            <option value="magma">Magma</option>
            <option value="jet">Jet</option>
            <option value="hot">Hot</option>
            <option value="cool">Cool</option>
        </select>
        <div id="colorbar"></div>
        <div class="label">{global_max:.2f}</div>
        <div style="height:160px;"></div>
        <div class="label">{global_min:.2f}</div>
    </div>

    <div id="controls" style="display:none;">
        <div id="frameInfo">Frame 1 / {len(frames_to_export)}</div>
        <div>
            <button onclick="firstFrame()">First</button>
            <button onclick="prevFrame()">Previous</button>
            <button id="playBtn" onclick="togglePlay()">Play</button>
            <button onclick="nextFrame()">Next</button>
            <button onclick="lastFrame()">Last</button>
            <button onclick="recenterModel()">Recenter</button>
            <button id="palmModeBtn" onclick="togglePalmMode()">Palm Mode: ON</button>
        </div>
        <input type="range" id="frameSlider" min="0" max="{len(frames_to_export)-1}" value="0" oninput="onSliderChange()">
        <div style="font-size:12px; margin-top:10px;">Desktop preview: drag to rotate, scroll to zoom</div>
    </div>

    <div id="arNotice" style="display:none;">
        Enable Hand Tracking in Quest Settings &gt; Movement Tracking, then tap "Enter AR" below
    </div>

    <script type="importmap">
    {{"imports": {{"three": "https://unpkg.com/three@0.160.0/build/three.module.js",
    "three/addons/": "https://unpkg.com/three@0.160.0/examples/jsm/"}}}}
    </script>

    <script type="module">
        import * as THREE from 'three';
        import {{ OrbitControls }} from 'three/addons/controls/OrbitControls.js';
        import {{ ARButton }} from 'three/addons/webxr/ARButton.js';

        const allMeta = {meta_json};
        console.log('Metadata loaded, frames:', allMeta.frames.length, 'property:', allMeta.property_name);

        const COLORMAPS = {{
            viridis: [[0.267,0.005,0.329],[0.283,0.141,0.458],[0.254,0.265,0.530],[0.164,0.471,0.558],[0.135,0.659,0.518],[0.478,0.821,0.318],[0.993,0.906,0.144]],
            plasma: [[0.050,0.030,0.528],[0.417,0.001,0.658],[0.693,0.165,0.565],[0.881,0.393,0.383],[0.987,0.810,0.145]],
            inferno: [[0.001,0.000,0.014],[0.144,0.047,0.338],[0.502,0.092,0.472],[0.855,0.258,0.283],[0.988,0.553,0.111]],
            magma: [[0.001,0.000,0.014],[0.129,0.046,0.285],[0.429,0.101,0.478],[0.670,0.188,0.427],[0.892,0.364,0.279],[0.987,0.810,0.145]],
            jet: t => {{const r=Math.max(0,Math.min(1,1.5-Math.abs(4*t-3))),g=Math.max(0,Math.min(1,1.5-Math.abs(4*t-2))),b=Math.max(0,Math.min(1,1.5-Math.abs(4*t-1)));return[r,g,b]}},
            hot: t => [Math.min(1,3*t),Math.min(1,Math.max(0,3*t-1)),Math.min(1,Math.max(0,3*t-2))],
            cool: t => [t,1-t,1]
        }};

        function getColor(cmap,t){{if(typeof cmap==='function')return cmap(t);const i=Math.min(Math.floor(t*(cmap.length-1)),cmap.length-2),f=t*(cmap.length-1)-i;return cmap[i].map((v,k)=>v+(cmap[i+1][k]-v)*f)}}

        let scene,camera,renderer,controls,particleSystem,modelGroup;
        let currentCmap='{default_cmap}',currentFrameIdx=0,isPlaying=false,playInterval,isLoadingFrame=false;

        const DESIRED_WORLD_SIZE = 0.35;
        const DEFAULT_POSITION = new THREE.Vector3(0, 1.2, -0.45);
        const initialScale = DESIRED_WORLD_SIZE / allMeta.size;
        const baseHalfExtent = allMeta.size * 0.5;

        scene = new THREE.Scene();
        scene.background = null;

        camera = new THREE.PerspectiveCamera(60, window.innerWidth/window.innerHeight, 0.01, 100);

        renderer = new THREE.WebGLRenderer({{ antialias:true, alpha:true }});
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setClearColor(0x000000, 0);
        renderer.xr.enabled = true;
        renderer.xr.setReferenceSpaceType('local-floor');
        document.getElementById('container').appendChild(renderer.domElement);

        document.body.appendChild(ARButton.createButton(renderer, {{
            optionalFeatures: ['hand-tracking', 'local-floor', 'dom-overlay'],
            domOverlay: {{ root: document.body }}
        }}));

        controls = new OrbitControls(camera, renderer.domElement);
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        camera.position.set(0, 1.2, 0.6);
        controls.target.copy(DEFAULT_POSITION);
        controls.update();

        scene.add(new THREE.AmbientLight(0xffffff, 1.0));
        const light = new THREE.DirectionalLight(0xffffff, 0.6);
        light.position.set(1, 2, 1);
        scene.add(light);

        modelGroup = new THREE.Group();
        modelGroup.position.copy(DEFAULT_POSITION);
        modelGroup.scale.setScalar(initialScale);
        scene.add(modelGroup);

        const refGrid = new THREE.GridHelper(allMeta.size * 1.2, 10, 0x0066ff, 0x003366);
        refGrid.position.y = -allMeta.size * 0.5;
        modelGroup.add(refGrid);
        modelGroup.add(new THREE.AxesHelper(allMeta.size * 0.3));

        async function updateFrame(){{
            if (isLoadingFrame) return;
            isLoadingFrame = true;
            const meta = allMeta.frames[currentFrameIdx];
            document.getElementById('particleInfo').innerHTML = `Loading ${{meta.count.toLocaleString()}} particles...`;

            const buf = await (await fetch(meta.file)).arrayBuffer();
            const data = new Float32Array(buf);
            const count = meta.count;

            if (particleSystem) modelGroup.remove(particleSystem);

            const geo = new THREE.BufferGeometry();
            const pos = new Float32Array(count*3);
            const col = new Float32Array(count*3);
            const cmap = COLORMAPS[currentCmap];

            for (let i=0; i<count; i++) {{
                pos[i*3]   = data[i*4]   - allMeta.center[0];
                pos[i*3+1] = data[i*4+1] - allMeta.center[1];
                pos[i*3+2] = data[i*4+2] - allMeta.center[2];
                const c = getColor(cmap, data[i*4+3]);
                col[i*3]=c[0]; col[i*3+1]=c[1]; col[i*3+2]=c[2];
            }}

            geo.setAttribute('position', new THREE.BufferAttribute(pos,3));
            geo.setAttribute('color', new THREE.BufferAttribute(col,3));

            const mat = new THREE.PointsMaterial({{
                size: allMeta.size * 0.008 * modelGroup.scale.x,
                vertexColors: true,
                sizeAttenuation: true
            }});
            particleSystem = new THREE.Points(geo, mat);
            modelGroup.add(particleSystem);

            document.getElementById('particleInfo').innerHTML =
                `Particles: ${{count.toLocaleString()}}<br>Property: {PROPERTY_NAME}<br>Range: ${{meta.property_min.toFixed(2)}} - ${{meta.property_max.toFixed(2)}}`;
            document.getElementById('frameInfo').textContent =
                `Frame ${{currentFrameIdx+1}} / ${{allMeta.frames.length}} (Step ${{meta.frame}})`;
            document.getElementById('frameSlider').value = currentFrameIdx;

            isLoadingFrame = false;
        }}

        function updateGradient(){{
            const cmap = COLORMAPS[currentCmap]; let grad = 'linear-gradient(to bottom';
            for (let i=0; i<=100; i++) {{ const c=getColor(cmap,1-i/100); grad+=`,rgb(${{Math.round(c[0]*255)}},${{Math.round(c[1]*255)}},${{Math.round(c[2]*255)}}) ${{i}}%` }}
            grad += ')'; document.getElementById('colorbar').style.background = grad;
        }}

        window.updateColormap = () => {{ currentCmap = document.getElementById('cmap').value; updateFrame(); updateGradient(); }};
        window.firstFrame = () => {{ currentFrameIdx = 0; updateFrame(); }};
        window.lastFrame = () => {{ currentFrameIdx = allMeta.frames.length-1; updateFrame(); }};
        window.nextFrame = () => {{ currentFrameIdx = (currentFrameIdx+1)%allMeta.frames.length; updateFrame(); }};
        window.prevFrame = () => {{ currentFrameIdx = (currentFrameIdx-1+allMeta.frames.length)%allMeta.frames.length; updateFrame(); }};
        window.onSliderChange = () => {{ currentFrameIdx = parseInt(document.getElementById('frameSlider').value); updateFrame(); }};
        window.recenterModel = () => {{
            modelGroup.position.copy(DEFAULT_POSITION);
            modelGroup.scale.setScalar(initialScale);
            modelGroup.rotation.set(0,0,0);
            if (particleSystem) particleSystem.material.size = allMeta.size*0.008*modelGroup.scale.x;
        }};
        window.togglePlay = () => {{
            if (isPlaying) {{
                clearInterval(playInterval);
                document.getElementById('playBtn').textContent='Play';
                isPlaying=false;
            }} else {{
                playInterval = setInterval(window.nextFrame, 700);
                document.getElementById('playBtn').textContent='Pause';
                isPlaying=true;
            }}
        }};

        document.addEventListener('keydown', e => {{
            if (e.key==='ArrowLeft') window.prevFrame();
            if (e.key==='ArrowRight') window.nextFrame();
            if (e.key===' ') {{ e.preventDefault(); window.togglePlay(); }}
        }});

        const hand1 = renderer.xr.getHand(0);
        const hand2 = renderer.xr.getHand(1);
        scene.add(hand1);
        scene.add(hand2);

        function makePinchMarker(){{
            const m = new THREE.Mesh(
                new THREE.SphereGeometry(0.008, 12, 12),
                new THREE.MeshBasicMaterial({{ color: 0x00ffff, transparent:true, opacity:0.85 }})
            );
            m.visible = false;
            scene.add(m);
            return m;
        }}
        const marker1 = makePinchMarker();
        const marker2 = makePinchMarker();

        const PINCH_THRESHOLD = 0.028;
        const _tmpThumb = new THREE.Vector3();
        const _tmpIndex = new THREE.Vector3();

        function getPinchPoint(hand, out){{
            const thumb = hand.joints['thumb-tip'];
            const index = hand.joints['index-finger-tip'];
            if (!thumb || !index) return null;
            thumb.getWorldPosition(_tmpThumb);
            index.getWorldPosition(_tmpIndex);
            if (_tmpThumb.distanceTo(_tmpIndex) > PINCH_THRESHOLD) return null;
            out.addVectors(_tmpThumb, _tmpIndex).multiplyScalar(0.5);
            return out;
        }}

        const pinch1 = new THREE.Vector3();
        const pinch2 = new THREE.Vector3();
        let singleGrab = null;
        let twoHandGrab = null;
        const worldPos = new THREE.Vector3();

        let hand1Handedness = null, hand2Handedness = null;
        hand1.addEventListener('connected', e => {{ hand1Handedness = e.data.handedness; }});
        hand1.addEventListener('disconnected', () => {{ hand1Handedness = null; }});
        hand2.addEventListener('connected', e => {{ hand2Handedness = e.data.handedness; }});
        hand2.addEventListener('disconnected', () => {{ hand2Handedness = null; }});

        let palmModeEnabled = true;
        window.togglePalmMode = () => {{
            palmModeEnabled = !palmModeEnabled;
            document.getElementById('palmModeBtn').textContent = `Palm Mode: ${{palmModeEnabled ? 'ON' : 'OFF'}}`;
        }};

        const WORLD_UP = new THREE.Vector3(0,1,0);
        const PALM_UP_THRESHOLD = 0.55;
        const PALM_HOVER_OFFSET = 0.05;
        const PALM_FOLLOW_LERP = 0.18;
        const PALM_NORMAL_SIGN = {{ left: 1, right: -1 }};

        const _palmWrist = new THREE.Vector3(), _palmIndex = new THREE.Vector3(), _palmPinky = new THREE.Vector3();
        const _palmPos = new THREE.Vector3(), _palmNormal = new THREE.Vector3();
        const _palmDirA = new THREE.Vector3(), _palmDirB = new THREE.Vector3();

        function getPalmPose(hand, handedness){{
            const wrist = hand.joints['wrist'];
            const indexMeta = hand.joints['index-finger-metacarpal'];
            const pinkyMeta = hand.joints['pinky-finger-metacarpal'];
            if (!wrist || !indexMeta || !pinkyMeta || !handedness) return null;

            wrist.getWorldPosition(_palmWrist);
            indexMeta.getWorldPosition(_palmIndex);
            pinkyMeta.getWorldPosition(_palmPinky);

            _palmPos.addVectors(_palmIndex, _palmPinky).multiplyScalar(0.5).add(_palmWrist).multiplyScalar(0.5);

            _palmDirA.subVectors(_palmIndex, _palmWrist);
            _palmDirB.subVectors(_palmPinky, _palmWrist);
            _palmNormal.crossVectors(_palmDirA, _palmDirB).normalize();
            _palmNormal.multiplyScalar(PALM_NORMAL_SIGN[handedness] || 1);

            return {{ position: _palmPos, normal: _palmNormal }};
        }}

        function updatePalmFollow(){{
            if (!palmModeEnabled || singleGrab || twoHandGrab) return false;
            for (const [hand, handedness] of [[hand1,hand1Handedness],[hand2,hand2Handedness]]) {{
                const pose = getPalmPose(hand, handedness);
                if (pose && pose.normal.dot(WORLD_UP) > PALM_UP_THRESHOLD) {{
                    const target = pose.position.clone().addScaledVector(pose.normal, PALM_HOVER_OFFSET);
                    modelGroup.position.lerp(target, PALM_FOLLOW_LERP);
                    return true;
                }}
            }}
            return false;
        }}

        function updateHandInteraction(){{
            if (!renderer.xr.isPresenting) return;

            const onPalm = updatePalmFollow();

            const p1 = getPinchPoint(hand1, pinch1);
            const p2 = getPinchPoint(hand2, pinch2);
            marker1.visible = !!p1; if (p1) marker1.position.copy(p1);
            marker2.visible = !!p2; if (p2) marker2.position.copy(p2);

            modelGroup.getWorldPosition(worldPos);
            const worldHalfExtent = baseHalfExtent * modelGroup.scale.x;

            if (onPalm) return;

            if (p1 && p2) {{
                const dist = p1.distanceTo(p2);
                const mid = new THREE.Vector3().addVectors(p1,p2).multiplyScalar(0.5);
                const angle = Math.atan2(p2.z - p1.z, p2.x - p1.x);
                const nearModel = Math.min(p1.distanceTo(worldPos), p2.distanceTo(worldPos)) < worldHalfExtent * 1.6;

                if (!twoHandGrab && nearModel) {{
                    twoHandGrab = {{
                        startDist: dist, startAngle: angle,
                        startScale: modelGroup.scale.x, startRotY: modelGroup.rotation.y,
                        startMid: mid.clone(), startPos: modelGroup.position.clone()
                    }};
                }} else if (twoHandGrab) {{
                    const newScale = THREE.MathUtils.clamp(twoHandGrab.startScale * (dist / twoHandGrab.startDist), initialScale*0.15, initialScale*6);
                    modelGroup.scale.setScalar(newScale);
                    modelGroup.rotation.y = twoHandGrab.startRotY + (angle - twoHandGrab.startAngle);
                    modelGroup.position.copy(twoHandGrab.startPos).add(new THREE.Vector3().subVectors(mid, twoHandGrab.startMid));
                    if (particleSystem) particleSystem.material.size = allMeta.size*0.008*modelGroup.scale.x;
                }}
                singleGrab = null;
            }} else {{
                twoHandGrab = null;
                const active = p1 ? {{ which:1, pt:p1 }} : (p2 ? {{ which:2, pt:p2 }} : null);

                if (active) {{
                    if (!singleGrab && active.pt.distanceTo(worldPos) < worldHalfExtent * 1.4) {{
                        singleGrab = {{ which: active.which, offset: new THREE.Vector3().subVectors(modelGroup.position, active.pt) }};
                    }} else if (singleGrab && singleGrab.which === active.which) {{
                        modelGroup.position.copy(active.pt).add(singleGrab.offset);
                    }}
                }} else {{
                    singleGrab = null;
                }}
            }}
        }}

        renderer.xr.addEventListener('sessionstart', () => {{ controls.enabled = false; }});
        renderer.xr.addEventListener('sessionend', () => {{ controls.enabled = true; singleGrab=null; twoHandGrab=null; }});

        async function init(){{
            document.getElementById('cmap').value = currentCmap;
            await updateFrame();
            updateGradient();
            document.getElementById('loading').style.display = 'none';
            document.getElementById('info').style.display = 'block';
            document.getElementById('colormap').style.display = 'block';
            document.getElementById('controls').style.display = 'block';
            document.getElementById('arNotice').style.display = 'block';
        }}
        init();

        renderer.setAnimationLoop(() => {{
            updateHandInteraction();
            if (controls.enabled) controls.update();
            renderer.render(scene, camera);
        }});
        window.addEventListener('resize', () => {{
            camera.aspect = window.innerWidth/window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }});

        console.log('READY - AR + hand tracking, streaming frames on demand');
    </script>
</body>
</html>"""

html_file = os.path.join(output_dir, f'{BASENAME}_ar.html')
with open(html_file, 'w', encoding='utf-8') as f:
    f.write(html)

print(f"\n{'='*70}")
print(f"AR HAND-TRACKING VIEWER CREATED!")
print(f"{'='*70}")
print(f"\nSaved: {html_file}")
print(f"Plus {len(frames_meta)} .bin files in the same folder ({bin_total_bytes/1024/1024:.1f} MB total)")
print(f"\nIMPORTANT: host the WHOLE '{os.path.basename(output_dir)}' folder together")
print(f"-- the .html file needs its .bin files sitting right next to it.")
print(f"\nColored by: {PROPERTY_NAME}")
print(f"Next step: upload that folder over HTTPS and open the .html link")
print(f"in the Meta Quest Browser, then tap 'Enter AR'.")
print(f"{'='*70}")