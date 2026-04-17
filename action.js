// action.js: action-minimizing solver for choreographic ring orbits.
//
// A "choreography" is a periodic orbit where all N rings trace the same
// closed path q(τ) in the (y,z) cross-section with equal time offsets
// τ + k/N for k = 0..N-1.  This file handles ONE orbit at a time; for
// systems with multiple distinct orbits (e.g. dense.html's shells), create
// one ActionMinimizer per orbit.
//
// Fourier parameterization (τ ∈ [0,1) is normalized time):
//   y(τ) = y0 + Σn  a[n-1]·cos(2πnτ) + b[n-1]·sin(2πnτ)
//   z(τ) =      Σn  c[n-1]·cos(2πnτ) + d[n-1]·sin(2πnτ)
//
// Action: S = ∫₀¹ [(y'²+z'²)/(2T) - T·V_eff(y,z)] dτ
// EOM in τ-space: y'' = T²·fy_eff  (from ÿ = fy_eff in physical time t = T·τ)
// Residual: Ry = y'' - T²·fy_eff   (→ 0 at the solution)
// Gradient: ∂S/∂aₙ = -(1/TM)·Σᵢ Ry·cos(2πnτᵢ)
// Update:   aₙ += (α/T)·(1/M)·Σᵢ Ry·cos(2πnτᵢ)
//
// T is kept fixed during step() — updating T during the step caused a
// positive-feedback collapse (T→0).  A future improvement would be to
// include T in the backtracking update, which test_gauge.js explored.
//
// Choreographic symmetry note: harmonics n = k·N (k≥1) produce nonzero
// total momentum when the N rings are sampled at equal phase offsets.
// They are zeroed in trajToParams() (initialization) and can be zeroed
// again via projectChoreography() after run().  They are NOT zeroed inside
// step() because for N=2 doing so causes the y-residual gradient to vanish
// for the kept (odd) modes, stalling convergence at rms ≈ 0.19.
//
// Requires ring.js (for ringForce).
//
// Usage:
//   var am = new ActionMinimizer({
//     N:     number of choreographic rings,
//     rm:    mass per ring,
//     l:     azimuthal angular momentum per ring,
//     fixed: [{y, z, mass}, ...]  // fixed background rings (e.g. central sun),
//     NF:    Fourier harmonics (default 8),
//     M:     quadrature points   (default 64),
//   });
//   am.trajToParams(traj, T_steps, inc);  // initialize from a raw trajectory
//   // or: am.params = { y0, T_phys, a[NF], b[NF], c[NF], d[NF] };
//   var rms = am.step(alpha);             // one gradient step
//   var rms = am.run({ nSteps, alpha, eps, reportEvery, onProgress });
//   am.projectChoreography();            // zero n=kN modes → zero total momentum
//   var samples = am.getSamples();        // N × {y, z, vy, vz}

(function() {
"use strict";

var ringForce = window.ringForce;  // captured from ring.js

// ---- Constructor ----

function ActionMinimizer(opts) {
    this.N      = opts.N;
    this.rm     = opts.rm;
    this.l      = opts.l;
    this.fixed  = opts.fixed || [];
    this.NF     = opts.NF || 8;
    this.M      = opts.M  || 64;
    this.params = null;
}

// ---- Fourier evaluation ----

// Evaluate orbit position and derivatives at normalized time τ.
// Returns {y, z, yp, zp, ypp, zpp} (position; first and second τ-derivatives).
ActionMinimizer.prototype.evalPt = function(tau) {
    var p = this.params, NF = this.NF;
    var y = p.y0, z = 0, yp = 0, zp = 0, ypp = 0, zpp = 0;
    for (var n = 1; n <= NF; n++) {
        var w = 2 * Math.PI * n, wt = w * tau;
        var cv = Math.cos(wt), sv = Math.sin(wt);
        var an = p.a[n-1], bn = p.b[n-1], cn = p.c[n-1], dn = p.d[n-1];
        y    += an*cv + bn*sv;
        z    += cn*cv + dn*sv;
        yp   += w*(-an*sv + bn*cv);
        zp   += w*(-cn*sv + dn*cv);
        ypp  -= w*w*(an*cv + bn*sv);
        zpp  -= w*w*(cn*cv + dn*sv);
    }
    return { y:y, z:z, yp:yp, zp:zp, ypp:ypp, zpp:zpp };
};

// ---- Force evaluation ----

// Effective force at (y, z) at normalized time tau:
//   centrifugal (l²/y³) + fixed-ring forces + peer ring forces.
ActionMinimizer.prototype.forceAt = function(y, z, tau) {
    var l = this.l, N = this.N, rm = this.rm;
    if (y < 1e-6) y = 1e-6;   // guard against centrifugal singularity
    var fy = l*l / (y*y*y);
    var fz = 0.0;
    for (var fi = 0; fi < this.fixed.length; fi++) {
        var fx = this.fixed[fi];
        var f  = ringForce(fx.y, fx.z, y, z);
        fy += f.fy * fx.mass;
        fz += f.fz * fx.mass;
    }
    for (var k = 1; k < N; k++) {
        var pt2 = this.evalPt((tau + k/N) % 1);
        var f   = ringForce(pt2.y, pt2.z, y, z);
        fy += f.fy * rm;
        fz += f.fz * rm;
    }
    return { fy:fy, fz:fz };
};

// ---- Gradient step ----

// One Newton-preconditioned step on the action residuals with backtracking line search.
// Returns RMS |residual| / T²  (→ 0 at convergence).
//
// Preconditioner: diagonal of the residual Jacobian dG/dx:
//   dGₙ/daₙ ≈ (1/M) Σ [−(2πn)² + T²·3l²/y⁴] cos²(2πnτ)  (y-modes)
//   dGₙ/dcₙ ≈ (1/M) Σ [−(2πn)²]              cos²(2πnτ)  (z-modes)
// When centrifugal dominates (hAn > wn²/10, low n), Newton subtracts.
// Otherwise (near-cancellation or kinetic-dominated modes), Sobolev adds.
//
// After computing the update direction, a backtracking line search halves alpha
// up to 8 times until rms strictly decreases.  This prevents divergence regardless
// of the initial alpha — the gradient direction is always correct for small enough steps.
ActionMinimizer.prototype.step = function(alpha) {
    var p = this.params, NF = this.NF, M = this.M;
    var l = this.l;
    var T = p.T_phys, T2 = T * T;
    var gA = new Array(NF).fill(0), gB = new Array(NF).fill(0);
    var gC = new Array(NF).fill(0), gD = new Array(NF).fill(0);
    var hA = new Array(NF).fill(0), hB = new Array(NF).fill(0);
    var rmsRes = 0;

    for (var i = 0; i < M; i++) {
        var tau = i / M;
        var pt  = this.evalPt(tau);
        var f   = this.forceAt(pt.y, pt.z, tau);
        var Ry  = pt.ypp - T2 * f.fy;
        var Rz  = pt.zpp - T2 * f.fz;
        rmsRes += Ry*Ry + Rz*Rz;
        var cent = T2 * 3 * l * l / (pt.y * pt.y * pt.y * pt.y);
        for (var n = 1; n <= NF; n++) {
            var wt = 2*Math.PI*n*tau, cv = Math.cos(wt), sv = Math.sin(wt);
            var wn2 = (2*Math.PI*n) * (2*Math.PI*n);
            gA[n-1] += Ry*cv/M;  gB[n-1] += Ry*sv/M;
            gC[n-1] += Rz*cv/M;  gD[n-1] += Rz*sv/M;
            hA[n-1] += (cent - wn2) * cv*cv / M;
            hB[n-1] += (cent - wn2) * sv*sv / M;
        }
    }
    var rms0 = Math.sqrt(rmsRes / M) / T2;

    // Compute preconditioned update directions (scale independent of alpha).
    var dA = new Array(NF), dB = new Array(NF), dC = new Array(NF), dD = new Array(NF);
    for (var n = 1; n <= NF; n++) {
        var wn2 = (2*Math.PI*n) * (2*Math.PI*n);
        var hAn = hA[n-1], hBn = hB[n-1];
        dA[n-1] = ((hAn > wn2 * 0.1) ? -1/hAn : 2/wn2) * gA[n-1];
        dB[n-1] = ((hBn > wn2 * 0.1) ? -1/hBn : 2/wn2) * gB[n-1];
        dC[n-1] = (2/wn2) * gC[n-1];
        dD[n-1] = (2/wn2) * gD[n-1];
    }
    // NOTE: Zeroing harmonics n = k·N here would enforce zero-total-momentum by
    // construction (see getSamples / trajToParams), but for N=2 it zeroes the
    // gradient of the y-residual for odd modes (they're orthogonal to the
    // constrained subspace), causing the solver to stall at a high residual.
    // The constraint is instead applied as a post-processing step in the caller
    // after the AM has converged.  See getSamples() comment.

    // Save params for backtracking.
    var savedA = p.a.slice(), savedB = p.b.slice();
    var savedC = p.c.slice(), savedD = p.d.slice();

    // Backtracking: halve alpha up to 8 times until rms decreases.
    var a = alpha;
    for (var k = 0; k <= 8; k++) {
        for (var n = 0; n < NF; n++) {
            p.a[n] = savedA[n] + a*dA[n];  p.b[n] = savedB[n] + a*dB[n];
            p.c[n] = savedC[n] + a*dC[n];  p.d[n] = savedD[n] + a*dD[n];
        }
        var rms1 = 0;
        for (var i = 0; i < M; i++) {
            var tau = i / M, pt = this.evalPt(tau), f = this.forceAt(pt.y, pt.z, tau);
            var Ry = pt.ypp - T2*f.fy, Rz = pt.zpp - T2*f.fz;
            rms1 += Ry*Ry + Rz*Rz;
        }
        rms1 = Math.sqrt(rms1/M) / T2;
        if (rms1 < rms0) return rms1;   // improvement found, keep update
        if (k < 8) a *= 0.5;
    }
    // No improvement after 8 halvings — restore and return original rms.
    p.a = savedA;  p.b = savedB;  p.c = savedC;  p.d = savedD;
    return rms0;
};

// ---- Run loop ----

// Run up to nSteps iterations with step size alpha.
// Stops early if rms < eps.
// onProgress(iter, rms, T_phys) is called every reportEvery steps (and on stop).
// Returns final rms.
ActionMinimizer.prototype.run = function(opts) {
    var nSteps      = opts.nSteps      || 2000;
    var alpha       = opts.alpha       || 1.0;
    var eps         = opts.eps         || 1e-10;
    var reportEvery = opts.reportEvery || 0;
    var onProgress  = opts.onProgress  || null;
    var rms = 0;
    for (var i = 0; i <= nSteps; i++) {
        rms = this.step(alpha);
        if (reportEvery > 0 && onProgress && i % reportEvery === 0)
            onProgress(i, rms, this.params.T_phys);
        if (rms < eps) {
            if (onProgress) onProgress(i, rms, this.params.T_phys);
            break;
        }
    }
    return rms;
};

// ---- Choreographic symmetry projection ----

// For an N-ring equally-spaced choreography the Fourier modes at harmonics
// n = N, 2N, 3N, … are the *only* ones that produce a nonzero total momentum
// (because the N phase-shifted copies add constructively instead of cancelling).
// Zeroing them enforces:
//   Σ_k vy(k/N) = 0,  Σ_k vz(k/N) = 0   (zero total linear momentum)
//   Σ_k y(k/N)  = N·y₀,  Σ_k z(k/N)  = 0 (COM at orbit centre)
//
// This projection is NOT applied inside step() because for N=2 it sets the
// gradient of the y-residual to zero for the odd modes that are kept, stalling
// convergence.  Instead call this once after run() has converged, then call
// getSamples().
ActionMinimizer.prototype.projectChoreography = function() {
    var p = this.params, NF = this.NF, N = this.N;
    if (N < 2) return;
    for (var n = N; n <= NF; n += N) {
        p.a[n-1] = 0;  p.b[n-1] = 0;
        p.c[n-1] = 0;  p.d[n-1] = 0;
    }
};

// ---- Output ----

// Return N initial-condition samples for Cosmos2D / ring():
//   [{y, z, vy, vz}, ...]  in physical units.
// For zero total momentum, call projectChoreography() first.
ActionMinimizer.prototype.getSamples = function() {
    var T = this.params.T_phys, N = this.N, out = [];
    for (var k = 0; k < N; k++) {
        var pt = this.evalPt(k / N);
        out.push({ y:pt.y, z:pt.z, vy:pt.yp/T, vz:pt.zp/T });
    }
    return out;
};

// ---- Initialization ----

// Build Fourier params from a raw trajectory.
//   traj[i] = [y, z, vy_phys, vz_phys]  (recorded at equal multistep intervals)
//   T_steps: period in multistep units
//   inc:     physical time per multistep  (= increment / work)
// Sets this.params and returns it.
ActionMinimizer.prototype.trajToParams = function(traj, T_steps, inc) {
    var NF = this.NF, M = this.M;
    var ys = [], zs = [];
    for (var i = 0; i < M; i++) {
        var t = i*T_steps/M, ii = Math.floor(t), f = t - ii;
        var p = traj[ii], q = traj[ii+1] || traj[ii];
        ys.push(p[0] + f*(q[0]-p[0]));
        zs.push(p[1] + f*(q[1]-p[1]));
    }
    var y0 = 0;
    for (var i = 0; i < M; i++) y0 += ys[i];
    y0 /= M;
    var a = [], b = [], c = [], d = [];
    for (var n = 1; n <= NF; n++) {
        var an = 0, bn = 0, cn = 0, dn = 0;
        for (var i = 0; i < M; i++) {
            var cv = Math.cos(2*Math.PI*n*i/M), sv = Math.sin(2*Math.PI*n*i/M);
            an += (ys[i]-y0)*cv;  bn += (ys[i]-y0)*sv;
            cn += zs[i]*cv;       dn += zs[i]*sv;
        }
        a.push(2*an/M);  b.push(2*bn/M);
        c.push(2*cn/M);  d.push(2*dn/M);
    }
    // Phase-align: rotate τ so that b[0] ≈ 0 (y-oscillation peaks at τ=0).
    // With b₁=0 the orbit is symmetric about τ=0, making the centrifugal
    // Hessian corrections well-posed (no near-cancellation between cos²/sin²
    // weighted averages) and eliminating the factor-of-4 asymmetry between
    // cosine and sine modes that causes slow or divergent convergence.
    var tau0 = Math.atan2(b[0], a[0]) / (2 * Math.PI);   // rotate to align a₁
    if (Math.abs(b[0]) > 0.01 * Math.abs(a[0])) {
        var newA = [], newB = [], newC = [], newD = [];
        for (var n = 1; n <= NF; n++) {
            var phase = 2 * Math.PI * n * tau0;
            var cp = Math.cos(phase), sp = Math.sin(phase);
            newA.push( a[n-1]*cp + b[n-1]*sp);
            newB.push(-a[n-1]*sp + b[n-1]*cp);
            newC.push( c[n-1]*cp + d[n-1]*sp);
            newD.push(-c[n-1]*sp + d[n-1]*cp);
        }
        a = newA;  b = newB;  c = newC;  d = newD;
    }
    // Zero choreographic-symmetry-breaking harmonics (same constraint as in step()).
    if (this.N >= 2) {
        for (var n = this.N; n <= NF; n += this.N) {
            a[n-1] = 0;  b[n-1] = 0;
            c[n-1] = 0;  d[n-1] = 0;
        }
    }
    this.params = { y0:y0, T_phys:T_steps*inc, a:a, b:b, c:c, d:d };
    return this.params;
};

// ---- HTML output ----

// Serialize a ring() options object (including moons) to a complete HTML file
// that only needs ring.js — no action.js, no solver, no computation.
// ringOpts is the full options object passed to ring(), moons included.
function actionOutputHTML(title, ringOpts) {
    var lines = [];
    lines.push('ring(\'cv\', {');
    for (var k in ringOpts) {
        if (k === 'moons') continue;
        lines.push('    ' + k + ': ' + JSON.stringify(ringOpts[k]) + ',');
    }
    lines.push('    moons: [');
    var moons = ringOpts.moons;
    for (var i = 0; i < moons.length; i++) {
        var m = moons[i], props = [];
        for (var mk in m) {
            var v = m[mk];
            props.push(mk + ': ' + (typeof v === 'number' ? v : JSON.stringify(v)));
        }
        lines.push('        {' + props.join(', ') + '}' + (i < moons.length - 1 ? ',' : ''));
    }
    lines.push('    ],');
    lines.push('});');

    return [
        '<!DOCTYPE html>',
        '<html>',
        '<head>',
        '  <title>' + title + '<\/title>',
        '  <script src="ring.js"><\/script>',
        '<\/head>',
        '<body bgcolor="#000000" text="#ffffff" style="font-family:monospace;">',
        '<p><canvas id="cv" width="600" height="400" style="border:1px solid #ffffff;"><\/canvas><\/p>',
        '<script>',
    ].concat(lines).concat([
        '<\/script>',
        '<\/body>',
        '<\/html>',
        '',
    ]).join('\n');
}

// Trigger a browser download of the generated HTML.
// filename defaults to <current-page-basename>_out.html (e.g. dyson.html → dyson_out.html).
function actionDownloadHTML(title, ringOpts, filename) {
    if (!filename) {
        var base = window.location.pathname.split('/').pop().replace(/\.html?$/i, '');
        filename = base + '_out.html';
    }
    var html = actionOutputHTML(title, ringOpts);
    var blob = new Blob([html], { type: 'text/html' });
    var a    = document.createElement('a');
    a.href     = URL.createObjectURL(blob);
    a.download = filename;
    a.click();
    URL.revokeObjectURL(a.href);
}

window.ActionMinimizer    = ActionMinimizer;
window.actionOutputHTML   = actionOutputHTML;
window.actionDownloadHTML = actionDownloadHTML;

})();
