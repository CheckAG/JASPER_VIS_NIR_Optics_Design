#!/usr/bin/env python3
"""Q4 -- FORWARD-optimize a REAL holographic recording for the simulated concave grating.

The three-panel "best case" (~0.9 nm at f/10) fits three FREE polynomial coefficients c1,c2,c3.
But a real grating is not made by dialing three independent numbers -- it is made by exposing the
blank to TWO coherent point sources C, D at a recording-laser wavelength lam0. Those five numbers
(rC, gamma, rD, delta, lam0) produce a SPECIFIC, coupled groove function

      nu(P) = (1/lam0) * (uC - uD) . g_hat(P)          [lines/mm]     (Noda-Namioka)

so the achievable c1,c2,c3 are locked to each other, not free. A free polynomial can therefore
promise a correction no two-point exposure can actually make.

This script answers the honest question: among the recordings a real holographic bench can make
(a discrete set of real laser lines, and geometry within reach), which one MINIMISES the spectral
resolution, and how close does that get to the free-polynomial "best case"?

  - If recording ~ polynomial  -> the best-case number is physically realizable (good).
  - If recording  > polynomial  -> the polynomial over-promised; the real target is the recording.

Everything is traced with ONE local tracer (below); the groove density is the only thing that
changes between CLASSICAL (constant), POLYNOMIAL (nu_of), and HOLOGRAPHIC (nu_holo). The polynomial
path is cross-checked against vls_raytrace.evaluate so this script is anchored to the main engine.

Simulated design: 230 l/mm, alpha=12.34 deg, arms 137.357/126.75 mm, R=133 mm (spherical substrate;
the toroid is the independent throughput knob and barely moves the dispersion-direction resolution).
Operating point f/10 -- the Panel-C best-case point. Pure NumPy + SciPy. No Optiland.
"""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar
import vls_raytrace as R

# ---------------- simulated design + operating point ----------------
N0 = 230.0; ALPHA = 12.34; rA = 137.357; rB0 = 126.75; Rsph = 133.0
BAND = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]; lam_nm = np.array(BAND) * 1e3
SLIT = 20.0; PIX = 8.0
FNUM = 10.0; EPD = rA / FNUM                       # f/10 -> the Panel-C best-case aperture
a = np.radians(ALPHA)
A = np.array([0.0, rA * np.sin(a), rA * np.cos(a)])  # slit point
Oc = np.array([0.0, 0.0, Rsph])                      # sphere centre

# Real holographic recording-laser lines to try (nm): quadrupled/tripled YAG, HeCd, Kr, Ar-ion, diode
LASERS = [266.0, 325.0, 363.8, 413.1, 441.6, 488.0]
PLATE0 = 34.2                                      # ~constant plate scale (nm/mm), for the fast search objective

def src(r, ang_deg):
    t = np.radians(ang_deg); return np.array([0.0, r * np.sin(t), r * np.cos(t)])

def nu_holo(P, g, C, D, lam0_mm):                  # physical holographic groove density (l/mm)
    uC = P - C; uC = uC / np.linalg.norm(uC, axis=1, keepdims=True)
    uD = P - D; uD = uD / np.linalg.norm(uD, axis=1, keepdims=True)
    return np.sum((uC - uD) * g, axis=1) / lam0_mm

# ---------------- one local tracer; groove density is the only variable ----------------
def _rays(npup):
    rng = np.random.default_rng(7)                 # same seed/sampling convention as vls_raytrace
    th = rng.uniform(0, 2 * np.pi, npup); rr = np.sqrt(rng.uniform(0, 1, npup)) * EPD / 2
    xs = rr * np.cos(th); ys = rr * np.sin(th)
    Pt = np.stack([xs, ys, (xs * xs + ys * ys) / (2 * Rsph)], axis=1)
    d = Pt - A; d /= np.linalg.norm(d, axis=1, keepdims=True)
    e = A - Oc; ed = d @ e; disc = ed * ed - (e @ e - Rsph * Rsph)
    t = -ed + np.sqrt(np.clip(disc, 0, None)); P = A + t[:, None] * d; Nrm = (Oc - P) / Rsph
    yhat = np.array([0., 1., 0.])
    g = yhat - (np.sum(yhat * Nrm, axis=1, keepdims=True)) * Nrm; g /= np.linalg.norm(g, axis=1, keepdims=True)
    return P, Nrm, g, d

def detector(rB):
    des = R.Design("sim", N0, ALPHA, rA, rB0, Rsph, Rsph, EPD, BAND)
    return des.detector([0, 0, 0], rB)

def trace(lam, nu_fn, Csys, npup=2500):
    Cc, ndet, edisp, ecross = Csys; lam_mm = lam * 1e-3
    P, Nrm, g, d = _rays(npup)
    nu = nu_fn(P, g)
    din = d - (np.sum(d * Nrm, axis=1, keepdims=True)) * Nrm
    dt = din + (lam_mm * nu)[:, None] * g
    n2 = 1.0 - np.sum(dt * dt, axis=1); ok = n2 >= 0
    sign = -np.sign(np.sum(d * Nrm, axis=1))
    dout = dt + (sign * np.sqrt(np.clip(n2, 0, None)))[:, None] * Nrm
    denom = dout @ ndet; tt = ((Cc - P) @ ndet) / denom; H = P + tt[:, None] * dout
    U = (H - Cc) @ edisp; V = (H - Cc) @ ecross; m = ok & np.isfinite(U)
    return U[m], V[m]

def band_resolution(nu_fn, rB, npup=1500):
    Csys = detector(rB)
    cen = [np.mean(trace(l, nu_fn, Csys, 1200)[0]) for l in BAND]
    rows = []
    for i, l in enumerate(BAND):
        U, V = trace(l, nu_fn, Csys, npup); rms = (U - U.mean()).std() * 1000
        aa = max(0, i - 1); bb = min(5, i + 1)
        plate = abs((lam_nm[bb] - lam_nm[aa]) / (cen[bb] - cen[aa]))
        diff = 1.028 * l * (rB / EPD)
        dl = np.hypot(np.hypot(SLIT, 2.355 * rms), np.hypot(PIX, diff)) / 1000 * plate
        rows.append(dict(lam=int(l * 1e3), rms=rms, dl=dl))
    mean = float(np.mean([r['dl'] for r in rows])); worst = float(max(r['dl'] for r in rows))
    return mean, worst, rows

def best_focus(nu_fn):
    r = minimize_scalar(lambda rB: band_resolution(nu_fn, rB, 250)[0],
                        bounds=(0.9 * rB0, 1.12 * rB0), method="bounded")
    return r.x

def _fast_res(C, D, lam0, rB, npup=200):
    """Lean search objective: RMS at a given detector distance, fixed plate scale, no centroid pass."""
    Csys = detector(rB); out = []
    for l in BAND:
        U, _ = trace(l, lambda P, g: nu_holo(P, g, C, D, lam0), Csys, npup)
        if len(U) < 40: return 1e6
        rms = (U - U.mean()).std() * 1000
        diff = 1.028 * l * (rB / EPD)
        out.append(np.hypot(np.hypot(SLIT, 2.355 * rms), np.hypot(PIX, diff)) / 1000 * PLATE0)
    out = np.array(out); return out.mean() + 0.6 * out.max()

# ---------------- forward optimisation over a real recording ----------------
def optimise_recording(lam0_nm):
    """For a fixed laser line, search (gamma, rC, rD, rB); delta fixed by the pole-density constraint
    sin(delta) = sin(gamma) + N0*lam0. Detector distance rB is a search variable (no nested focus),
    so each objective call is a single band trace. Returns the best achievable resolution + geometry."""
    lam0 = lam0_nm * 1e-6; k = N0 * lam0
    def unpack(x):
        gamma, rC, rD, rB = x; s = np.sin(np.radians(gamma)) + k
        if abs(s) > 0.999: return None
        if not (60 < rC < 1200 and 60 < rD < 1200 and 0.9 * rB0 < rB < 1.12 * rB0): return None
        return gamma, rC, rD, np.degrees(np.arcsin(s)), rB
    def resfn(x):
        u = unpack(x)
        if u is None: return 1e6
        gamma, rC, rD, delta, rB = u
        return _fast_res(src(rC, gamma), src(rD, delta), lam0, rB)
    best = None
    for seed in ([5, 300, 300, rB0], [25, 250, 450, rB0], [-25, 450, 250, rB0]):
        r = minimize(resfn, seed, method="Nelder-Mead",
                     options=dict(xatol=1e-3, fatol=1e-2, maxiter=300))
        if best is None or r.fun < best.fun: best = r
    u = unpack(best.x)
    if u is None: return None
    gamma, rC, rD, delta, rBs = u
    C = src(rC, gamma); D = src(rD, delta)
    nu_fn = lambda P, g: nu_holo(P, g, C, D, lam0)
    rB = best_focus(nu_fn)                          # refine focus with the full centroid pass at the end
    m, w, rows = band_resolution(nu_fn, rB, 2000)
    gp = np.array([[0., 1., 0.]]); nup = nu_holo(np.zeros((1, 3)), gp, C, D, lam0)[0]
    return dict(lam0=lam0_nm, gamma=gamma, rC=rC, rD=rD, delta=delta, rB=rB,
                mean=m, worst=w, rows=rows, nup=nup, C=C, D=D, lam0_mm=lam0)

def main():
    print("=" * 78)
    print("  Q4: best REAL 2-point holographic recording for the simulated grating (f/10)")
    print("=" * 78)

    # --- reference: free-polynomial best case, from the main engine, same design + f/10 ---
    des = R.Design("sim", N0, ALPHA, rA, rB0, Rsph, Rsph, EPD, BAND)
    c = R.fit_vls(des)
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        poly = R.evaluate(des, c, "poly")
    print(f"\nFREE-POLYNOMIAL best case (fitted c1,c2,c3 = "
          f"[{c[0]:.4g}, {c[1]:.4g}, {c[2]:.4g}]):")
    print(f"  Dlam mean {poly['mean']:.2f}  worst {poly['worst']:.2f} nm   <- the number a real recording must reach")
    # cross-check: reproduce the polynomial number with THIS script's tracer
    nu_poly = lambda P, g: R.nu_of(P[:, 1], N0, *c)
    rBp = best_focus(nu_poly); mp, wp, _ = band_resolution(nu_poly, rBp, 2000)
    print(f"  (this script's tracer on the same c: mean {mp:.2f} worst {wp:.2f} nm -> tracer anchored to engine)")

    # --- forward-optimise over each real laser line ---
    print(f"\nScanning {len(LASERS)} real recording-laser lines...")
    results = []
    for lam0 in LASERS:
        s = optimise_recording(lam0)
        if s is None: continue
        results.append(s)
        print(f"  lam0={lam0:6.1f} nm  ->  C=(r{s['rC']:6.1f}mm,{s['gamma']:+6.2f}d) "
              f"D=(r{s['rD']:6.1f}mm,{s['delta']:+6.2f}d)  pole {s['nup']:.1f} l/mm  "
              f"Dlam mean {s['mean']:.2f} worst {s['worst']:.2f} nm")

    best = min(results, key=lambda s: s['mean'] + 0.6 * s['worst'])
    print("\n" + "-" * 78)
    print(f"BEST REAL RECORDING: lam0={best['lam0']:.1f} nm | "
          f"C=(r{best['rC']:.1f} mm, {best['gamma']:+.2f} deg)  D=(r{best['rD']:.1f} mm, {best['delta']:+.2f} deg)")
    print(f"  pole density {best['nup']:.1f} l/mm (target {N0}) | detector r'={best['rB']:.1f} mm")
    print(f"  achievable  Dlam mean {best['mean']:.2f}  worst {best['worst']:.2f} nm")
    print(f"  free-poly   Dlam mean {poly['mean']:.2f}  worst {poly['worst']:.2f} nm")
    gap = best['mean'] - poly['mean']
    verdict = ("REALIZABLE: a real 2-point exposure reaches the polynomial best case"
               if gap < 0.15 else
               "PARTIAL: a simple 2-point exposure falls short of the free polynomial by "
               f"{gap:.2f} nm -> needs aspheric/CGH recording for the rest")
    print(f"  verdict: {verdict}")

    import json
    json.dump(dict(design="simulated", fnum=FNUM,
                   poly=dict(c=list(c), mean=poly['mean'], worst=poly['worst']),
                   best_recording={k: best[k] for k in ('lam0','gamma','rC','rD','delta','rB','mean','worst','nup')},
                   scan=[{k: s[k] for k in ('lam0','gamma','rC','rD','delta','mean','worst','nup')} for s in results]),
              open("recording_optimize_sim.json", "w"), indent=1)
    make_figure(best, c, poly, results)
    spot_diagram(best)
    print("\nsaved jasper_recording_optimize.png, jasper_recording_spots.png and recording_optimize_sim.json")

def spot_diagram(best, fname="jasper_recording_spots.png"):
    """Per-wavelength ray-traced spots for the BEST REAL RECORDING — the physical holographic grooves
    nu_holo(C,D,lam0) applied directly (NOT a polynomial), at its own best-focus detector."""
    Csys = detector(best['rB'])
    nu_fn = lambda P, g: nu_holo(P, g, best['C'], best['D'], best['lam0_mm'])
    fig, axes = plt.subplots(1, len(BAND), figsize=(2.05 * len(BAND), 2.6), sharey=True)
    for ax, (l, row) in zip(axes, zip(BAND, best['rows'])):
        U, V = trace(l, nu_fn, Csys, 2500)
        u = (U - U.mean()) * 1000.0; v = (V - V.mean()) * 1000.0     # mm -> um, centroid-relative
        ax.scatter(u, v, s=1.5, color="#1aa5cf", alpha=0.5, edgecolors="none")
        ax.set_title(f"{int(l * 1e3)} nm\nΔλ {row['dl']:.2f} nm", fontsize=8.5)
        ax.axhline(0, color="0.85", lw=0.5); ax.axvline(0, color="0.85", lw=0.5)
        ax.set_xlim(-40, 40); ax.tick_params(labelsize=7)
        ax.set_xlabel("dispersion (µm)", fontsize=7.5)
    axes[0].set_ylabel("cross-dispersion (µm)", fontsize=8)
    fig.suptitle(f"Simulated grating at f/10 — ray-traced spots for the best REAL recording "
                 f"(λ₀={best['lam0']:.1f} nm, physical ν(y))\n"
                 f"pole {best['nup']:.0f} l/mm · Δλ mean {best['mean']:.2f} / worst {best['worst']:.2f} nm "
                 f"(resolution = dispersion/horizontal width; the vertical astigmatic line is on the "
                 f"spherical blank — the toroid collapses it)",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(fname, dpi=150, bbox_inches="tight"); plt.close(fig)

def make_figure(best, c, poly, results):
    C_HOLO, C_POLY, C_ACC = "#1aa5cf", "#e8912a", "#2f8f5b"
    ys = np.linspace(-EPD / 2, EPD / 2, 220)
    P = np.stack([np.zeros_like(ys), ys, ys * ys / (2 * Rsph)], axis=1)
    Nrm = (Oc - P) / Rsph; yhat = np.array([0., 1., 0.])
    g = yhat - (np.sum(yhat * Nrm, axis=1, keepdims=True)) * Nrm; g /= np.linalg.norm(g, axis=1, keepdims=True)
    nu_h = nu_holo(P, g, best['C'], best['D'], best['lam0_mm'])
    nu_p = R.nu_of(ys, N0, *c)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Q4 - a REAL 2-point holographic recording vs the free-polynomial best case (simulated grating, f/10)",
                 fontsize=12.5, fontweight="bold", y=1.0)

    # left: groove profiles
    ax1.plot(ys, nu_h, lw=2.4, color=C_HOLO,
             label=f"physical holographic ν(y)\n(real recording, λ₀={best['lam0']:.1f} nm)")
    ax1.plot(ys, nu_p, lw=1.8, ls="--", color=C_POLY, label="free-polynomial ν(y) (unconstrained best case)")
    ax1.axhline(N0, color="0.7", lw=0.8, ls=":"); ax1.text(-EPD / 2, N0 + 0.03, "pole 230 l/mm", fontsize=8, color="0.5")
    ax1.set_xlabel("position across aperture  y (mm)"); ax1.set_ylabel("local groove density  ν (l/mm)")
    ax1.set_title("Groove profile: physical vs free", fontsize=10.5)
    ax1.grid(alpha=0.3); ax1.legend(fontsize=8.5, loc="best")

    # right: resolution -- best real recording vs polynomial, per wavelength
    lams = [r['lam'] for r in poly['rows']]
    ax2.plot(lams, [r['dlam'] for r in poly['rows']], "o--", color=C_POLY, lw=1.8, label=f"free polynomial (mean {poly['mean']:.2f} nm)")
    ax2.plot(lams, [r['dl'] for r in best['rows']], "s-", color=C_HOLO, lw=2.2, label=f"best real recording (mean {best['mean']:.2f} nm)")
    ax2.axhline(1.0, color=C_ACC, lw=1.0, ls=":"); ax2.text(505, 1.02, "1 nm", fontsize=8, color=C_ACC)
    ax2.set_xlabel("wavelength (nm)"); ax2.set_ylabel("resolution  Δλ (nm)")
    ax2.set_title("Resolution across the band", fontsize=10.5)
    ax2.grid(alpha=0.3); ax2.legend(fontsize=8.5, loc="best")
    ax2.set_ylim(0, max(1.3, max(r['dl'] for r in best['rows']) * 1.15))

    fig.tight_layout()
    fig.savefig("jasper_recording_optimize.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
