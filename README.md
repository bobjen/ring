# ring

JavaScript gravitational ring simulator.  Models the dynamics of massive rings
sharing a common axis, displayed as a YZ cross-section.

Gravitational force computation ported from `ringforce.cpp` by Bob Jenkins
(2017, public domain).  JavaScript simulator written by Claude 4.5 (February
2026) with interface identical to `orbit.js`.

Documentation and live demos: https://burtleburtle.net/bob/js/ring.html

## How it works

Each ring is modeled as 8,192 equally-spaced point masses arranged in a circle
around the X=0, Y=0 axis.  The gravitational force between two rings is
computed analytically using Chebyshev interpolation of elliptic integrals,
avoiding the need to simulate all 8,192² pairs explicitly.

Integration uses the same explicit symmetric multistep method as `orbit.js`
(1–15 prior accelerations, default 8).

Angular momentum `l = Y² dφ/dt` is stored as a conserved quantity rather than
tracking azimuthal velocity directly, preventing cumulative errors from causing
drift.

## Coordinate system

```
Y = radial distance from axis  (always > 0, maps to screen X)
Z = axial position              (can be negative, maps to screen Y)
```

The X axis (azimuthal) is integrated out by the ring symmetry.

## Usage

```html
<script src="ring.js"></script>
<canvas id="sim" width="600" height="400"></canvas>
<script>
ring("sim", {
    background: "#000",
    scale: 50,
    increment: 0.06,
    work: 30,
    trail: true,
    moons: [
        // Central sun (fixed on axis)
        { y: 0.01, z: 0, mass: 100, color: "#ffff00", radius: 4, fixed: true },
        // Orbital ring with Keplerian angular momentum
        { y: 5, z: 0, mass: 1, color: "#00aaff", radius: 4,
          l: Math.sqrt(100 * 5) },
    ]
});
</script>
```

## Moon object fields

| Field | Default | Description |
|-------|---------|-------------|
| `y` | 0 | Radial distance from axis |
| `z` | 0 | Axial position |
| `vy` | 0 | Radial velocity |
| `vz` | 0 | Axial velocity |
| `l` | 0 | Conserved angular momentum (use instead of `vx`) |
| `vx` | — | Azimuthal speed (alternative to `l`; sets `l = y * vx`) |
| `mass` | 0 | Mass (0 = massless tracer) |
| `color` | `"#ffffff"` | Display color |
| `radius` | 3 | Display radius in pixels |
| `fixed` | false | If true, radial position is held constant |

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `moons` | (required) | Array of ring objects |
| `increment` | 1.0 | Simulated time per displayed frame |
| `work` | 20 | Integration steps per frame |
| `points` | 8 | Multistep order (1–15) |
| `framerate` | 50 | Milliseconds between frames |
| `scale` | 50 | Pixels per unit length |
| `ymargin` | 20 | Pixels from left edge to Y=0 |
| `background` | `"#000000"` | Canvas background color |
| `trail` | false | Draw trails |
| `traillen` | 200 | Maximum trail length in points |
| `fade` | — | Trail fade alpha per frame (e.g. 0.05) |
| `fadeto` | background | Color to fade toward |
| `stop` | false | Start paused (click canvas to toggle) |
| `lifetime` | 0 (forever) | Stop after this many simulated time units |

## Building

No build step — `ring.js` is a single self-contained script.  Open `test.html`
in a browser to see the included demos.
