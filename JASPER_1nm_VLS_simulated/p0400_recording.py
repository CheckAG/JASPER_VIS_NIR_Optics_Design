#!/usr/bin/env python3
"""lambda_0 back-solve for the Shimadzu P0400-01: which two-point holographic recording reproduces it?

Everything fixed from the datasheet (see validate_p0400.py): 400 l/mm at the pole, alpha = 20 deg,
r = 96.9 mm, R = 99.61 mm (blank sag), detector on the datasheet plane, 400-1100 nm, f/5.3.

For each candidate recording-laser line lambda_0 we search the recording geometry -- two point sources
C = (r_C, gamma) and D = (r_D, delta) -- with the pole density pinned by physics,
      sin(delta) = sin(gamma) + nu0 * lambda_0 ,
and trace the PHYSICAL groove function nu(P) = (1/lambda_0)(u_C - u_D).g(P) directly (no polynomial).
We then compare (a) the resolution each recording reaches with the best free polynomial (2.02 nm), and
(b) the groove profile it writes with the fitted polynomial nu = 400 + c1 y + c2 y^2 + c3 y^3.

Honest limit: optical performance constrains only the groove profile over the illuminated aperture.
Different (lambda_0, geometry) pairs can write nearly the same profile, so lambda_0 may not be uniquely
identifiable from datasheet data -- the scan below shows whether it is.
"""
import numpy as np, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, minimize
import validate_p0400 as V

N0 = 400.0; ALPHA = 20.0; rA = 96.9; RS = V.RSUB; EPD = rA / 5.3
BAND = V.BAND; CS = V.CS; lam_nm = np.array(BAND) * 1e3
a = np.radians(ALPHA); A = np.array([0.0, rA * np.sin(a), rA * np.cos(a)]); Oc = np.array([0.0, 0.0, RS])
LASERS = [266.0, 325.0, 351.1, 413.1, 441.6, 457.9, 488.0]
C_POLY = [0.2528, -0.01713, -3.0e-7]          # best free polynomial on the same plane (validate_models)

def src(r, ang): t = np.radians(ang); return np.array([0.0, r * np.sin(t), r * np.cos(t)])

def nu_holo(P, g, C, D, lam0):
    uC = P - C; uC /= np.linalg.norm(uC, axis=1, keepdims=True)
    uD = P - D; uD /= np.linalg.norm(uD, axis=1, keepdims=True)
    return np.sum((uC - uD) * g, axis=1) / lam0

_cache = {}
def rays(npup):
    if npup in _cache: return _cache[npup]
    rng = np.random.default_rng(7)
    th = rng.uniform(0, 2 * np.pi, npup); rr = np.sqrt(rng.uniform(0, 1, npup)) * EPD / 2
    xs, ys = rr * np.cos(th), rr * np.sin(th)
    Pt = np.stack([xs, ys, (xs * xs + ys * ys) / (2 * RS)], axis=1)
    d = Pt - A; d /= np.linalg.norm(d, axis=1, keepdims=True)
    e = A - Oc; ed = d @ e; t = -ed + np.sqrt(np.clip(ed * ed - (e @ e - RS * RS), 0, None))
    P = A + t[:, None] * d; Nrm = (Oc - P) / RS
    g = np.array([0., 1., 0.]) - Nrm[:, 1:2] * Nrm; g /= np.linalg.norm(g, axis=1, keepdims=True)
    _cache[npup] = (P, Nrm, g, d); return _cache[npup]

def trace(lam, nu_fn, npup):
    Cc, ndet, ed, ex = CS; P, Nrm, g, d = rays(npup)
    din = d - np.sum(d * Nrm, axis=1, keepdims=True) * Nrm
    dt = din + (lam * 1e-3 * nu_fn(P, g))[:, None] * g
    n2 = 1 - np.sum(dt * dt, axis=1); ok = n2 >= 0
    dout = dt + (-np.sign(np.sum(d * Nrm, axis=1)) * np.sqrt(np.clip(n2, 0, None)))[:, None] * Nrm
    H = P + (((Cc - P) @ ndet) / (dout @ ndet))[:, None] * dout
    U = (H - Cc) @ ed; V_ = (H - Cc) @ ex; m = ok & np.isfinite(U)
    return U[m], V_[m]

def resolution(nu_fn):
    cen = [np.mean(trace(l, nu_fn, 1200)[0]) for l in BAND]; out = []
    for i, l in enumerate(BAND):
        U, _ = trace(l, nu_fn, 2500); rms = (U - U.mean()).std() * 1e3
        i0, i1 = max(0, i - 1), min(len(BAND) - 1, i + 1)
        plate = abs((lam_nm[i1] - lam_nm[i0]) / (cen[i1] - cen[i0]))
        out.append(np.hypot(np.hypot(20, 2.355 * rms), np.hypot(8, 1.028 * l * V.plane_r(l) / EPD)) / 1e3 * plate)
    return np.array(out)

def geometry(x, lam0_nm):
    gamma, rC, rD = x; s = np.sin(np.radians(gamma)) + N0 * lam0_nm * 1e-6
    if abs(s) >= 0.99: return None
    return gamma, rC, rD, np.degrees(np.arcsin(s))

def cost(x, lam0_nm):
    geo = geometry(x, lam0_nm)
    if geo is None: return 1e6
    gamma, rC, rD, delta = geo; C, D = src(rC, gamma), src(rD, delta); lam0 = lam0_nm * 1e-6
    out = []
    for l in BAND:
        U, _ = trace(l, lambda P, g: nu_holo(P, g, C, D, lam0), 300)
        if len(U) < 50: return 1e6
        out.append((U - U.mean()).std() * 1e3)
    out = np.array(out); return out.mean() + 0.6 * out.max()

def solve(lam0_nm):
    bounds = [(-75, 75), (50, 2000), (50, 2000)]
    r = differential_evolution(cost, bounds, args=(lam0_nm,), seed=3, maxiter=70, popsize=18, tol=1e-6, polish=False)
    r = minimize(cost, r.x, args=(lam0_nm,), method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-4, maxiter=1500))
    gamma, rC, rD, delta = geometry(r.x, lam0_nm); C, D = src(rC, gamma), src(rD, delta); lam0 = lam0_nm * 1e-6
    dl = resolution(lambda P, g: nu_holo(P, g, C, D, lam0))
    # groove profile along the meridian and its cubic fit over the illuminated aperture
    y = np.linspace(-20, 20, 401); P = np.stack([np.zeros_like(y), y, RS - np.sqrt(RS * RS - y * y)], axis=1)
    Nrm = (Oc - P) / RS; g = np.array([0., 1., 0.]) - Nrm[:, 1:2] * Nrm; g /= np.linalg.norm(g, axis=1, keepdims=True)
    nu = nu_holo(P, g, C, D, lam0); m = np.abs(y) <= EPD / 2
    cf = np.polyfit(y[m], nu[m] - N0, 3)[::-1][1:]          # c1, c2, c3 (constant term ~0)
    return dict(lam0=lam0_nm, gamma=gamma, rC=rC, rD=rD, delta=delta, dl=dl.tolist(), mean=float(dl.mean()),
                worst=float(dl.max()), c=cf.tolist(), nu=nu.tolist(), y=y.tolist(), nu_pole=float(nu[200]))

def main():
    print("=" * 96); print("  lambda_0 back-solve: two-point holographic recording for Shimadzu P0400-01"); print("=" * 96)
    print(f"  fixed: 400 l/mm at pole, alpha 20 deg, r 96.9, R {RS:.2f} (blank), datasheet plane, 400-1100 nm, f/5.3")
    print(f"  target: best free polynomial c = {C_POLY} -> 2.02 nm mean\n")
    res = []
    for lam0 in LASERS:
        s = solve(lam0); res.append(s)
        print(f"  lambda0 {lam0:6.1f} nm | C=({s['rC']:7.1f} mm, {s['gamma']:+6.2f} deg)  D=({s['rD']:7.1f} mm, {s['delta']:+6.2f} deg)"
              f" | pole {s['nu_pole']:.2f} | c1 {s['c'][0]:+.4f} c2 {s['c'][1]:+.5f} | mean {s['mean']:.2f} worst {s['worst']:.2f} nm")
    json.dump(res, open("p0400_recording.json", "w"), indent=1)
    make_figure(res); print("\nsaved p0400_recording.png and p0400_recording.json")

def make_figure(res):
    cols = plt.cm.viridis(np.linspace(0.05, 0.9, len(res)))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 4.8))
    fig.suptitle("P0400-01 — two-point holographic recordings that reproduce the grating, per recording laser",
                 fontsize=12.5, fontweight="bold", y=1.02)
    y = np.array(res[0]["y"]); poly = N0 + C_POLY[0] * y + C_POLY[1] * y**2 + C_POLY[2] * y**3
    a1.axvspan(-EPD / 2, EPD / 2, color="0.93", zorder=0); a1.text(0, 386.5, "illuminated aperture (f/5.3)", ha="center", fontsize=8.5, color="0.4")
    a1.plot(y, poly, "k--", lw=2, label="best free polynomial", zorder=5)
    for s, c in zip(res, cols): a1.plot(y, s["nu"], color=c, lw=1.6, label=f"λ₀ {s['lam0']:.1f} nm")
    a1.set_xlabel("position across the 40 mm blank, y (mm)"); a1.set_ylabel("groove density ν (l/mm)")
    a1.set_title("Groove profile written by each recording", fontsize=10.5); a1.grid(alpha=0.3); a1.legend(fontsize=7.5, ncol=2)
    a1.set_ylim(385, 403)
    a2.bar(range(len(res)), [s["mean"] for s in res], color=cols)
    a2.axhline(2.02, color="k", ls="--", lw=1.2); a2.text(-0.45, 2.33, "dashed: best free polynomial, 2.02 nm", ha="left", fontsize=8.5)
    for i, s in enumerate(res): a2.text(i, s["mean"] + 0.04, f"{s['mean']:.2f}", ha="center", fontsize=8.5)
    a2.set_xticks(range(len(res))); a2.set_xticklabels([f"{s['lam0']:.1f}" for s in res])
    a2.set_xlabel("recording wavelength λ₀ (nm)"); a2.set_ylabel("mean resolution 400–1100 nm (nm)")
    a2.set_title("Resolution each recording reaches (f/5.3, datasheet plane)", fontsize=10.5); a2.grid(axis="y", alpha=0.3)
    a2.set_ylim(0, max(2.6, max(s["mean"] for s in res) * 1.15))
    fig.tight_layout(); fig.savefig("p0400_recording.png", dpi=150, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    main()
