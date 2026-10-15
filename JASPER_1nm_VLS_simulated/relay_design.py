#!/usr/bin/env python3
"""Fiber -> relay -> 25 um slit -> Shimadzu P0400-01 -> detector: an end-to-end, real-surface model.

Two catalog relays (prescriptions from the Edmund datasheets; glass + EFL/BFL verified in relay_optics.py):
  compact  #49-321 12.5 mm dia x 14 mm FL (N-BAF10/N-SF5, reversed) + #49-329 12.5 mm dia x 60 mm FL (N-BK7/N-SF5)
           magnification 60/14 = 4.29 -> NA 0.22 at the fiber becomes NA 0.051 at the slit (f/9.8)
  large    #49-352 25 mm dia x 30 mm FL (N-BAF10/N-SF10, reversed) + #49-361 25 mm dia x 125 mm FL (N-BK7/N-SF5)
           magnification 125/30 = 4.17 -> NA 0.053 at the slit (f/9.5)
Fiber: Edmund #23-692, 200 um core, NA 0.22 VIS-NIR, fully and uniformly filled (extended diffuse source).
Slit: Edmund #58-541, 25 +/-3 um x 3 mm. Grating: P0400-01 model from validate_p0400.py (f/10 solved chirp).
Detector: 8 x 200 um pixels on the datasheet plane. Recommended: no stop on the grating (the relay's NA
limits the beam); the f/10-stop case is reported for comparison.

Each ray is traced through the real lens surfaces, clipped by the slit, then traced through the concave VLS
grating to the detector. Per wavelength:
  relay+slit  fraction of the fiber's light that passes the slit and reaches the grating aperture
  on-pixel    fraction of that landing in the 200 um-tall pixel row
  system      product -- geometric only (coating, fiber, grating, detector efficiencies not included)
  dlam        package formula: slit (25 um), 2.355 x aberration RMS (traced from the slit centre under the
              relay's real illumination), 8 um pixel, diffraction; x plate scale on the datasheet plane
  dlam_fwhm   direct check: FWHM of the traced line profile (all slit rays on the pixel row, convolved
              with the 8 um pixel, diffraction added in quadrature) x plate scale
The slit sits at the axial position that maximises band-averaged system throughput.

Run:  python relay_design.py            # both relays (~4 min)
      python relay_design.py compact    # one relay
Outputs: relay_design_<name>.json, relay_performance_<name>.png, relay_spots_<name>.png, relay_layouts.png
"""
import sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import relay_optics as RO, vls_raytrace as R, validate_p0400 as V

CORE = 0.200; NAF = 0.22; SLIT_W = 0.025; SLIT_H = 3.0; NRAY = 300_000; PIXW = 8.0; PIXH = 200.0
BAND = V.BAND; LAM = np.array(BAND) * 1e3
CHIRP = json.load(open("validate_p0400.json"))["f/10"]["c"]
EPD10 = 96.9 / 10
COLS = plt.cm.turbo(np.linspace(0.08, 0.92, len(BAND)))
PAIRS = {"compact": dict(coll="49-321", foc="49-329", gap=5.0),
         "large":   dict(coll="49-352", foc="49-361", gap=10.0)}

def config(name):
    p = PAIRS[name]; C = dict(name=name, **p)
    C["COLL"] = RO.surfaces(p["coll"], reverse=True); C["FOC"] = RO.surfaces(p["foc"])
    C["ca_c"] = RO.LENS[p["coll"]]["ca"] / 2; C["ca_f"] = RO.LENS[p["foc"]]["ca"] / 2
    C["rad_c"] = RO.LENS[p["coll"]]["dia"] / 2; C["rad_f"] = RO.LENS[p["foc"]]["dia"] / 2
    C["D1"] = RO.paraxial(RO.surfaces(p["coll"]), 0.5876)[1]          # fiber at the collimator's back focus
    C["Z_FOC"] = C["D1"] + sum(t for _, t, _ in C["COLL"]) + p["gap"]
    C["Z_SLIT0"] = C["Z_FOC"] + sum(t for _, t, _ in C["FOC"]) + RO.paraxial(C["FOC"], 0.5876)[1]
    C["m"] = RO.LENS[p["foc"]]["efl"] / RO.LENS[p["coll"]]["efl"]
    C["label"] = (f"#{p['coll']} {RO.LENS[p['coll']]['efl']:.0f} mm + #{p['foc']} {RO.LENS[p['foc']]['efl']:.0f} mm "
                  f"({RO.LENS[p['coll']]['dia']:g} mm dia)")
    return C

def fiber_rays(N, seed=11):
    rng = np.random.default_rng(seed)
    r = CORE / 2 * np.sqrt(rng.uniform(0, 1, N)); th = rng.uniform(0, 2 * np.pi, N)
    k = NAF * np.sqrt(rng.uniform(0, 1, N)); ph = rng.uniform(0, 2 * np.pi, N)     # uniform radiance inside the NA
    return (np.stack([r * np.cos(th), r * np.sin(th), np.zeros(N)], 1),
            np.stack([k * np.cos(ph), k * np.sin(ph), np.sqrt(1 - k * k)], 1))

_rays = fiber_rays(NRAY); _cache = {}
def relay_out(C, lam):
    key = (C["name"], lam)
    if key not in _cache:
        P, D = _rays
        P, D, ok, _ = RO.trace(P, D, C["COLL"], C["D1"], lam, C["ca_c"])
        P, D, ok2, _ = RO.trace(P, D, C["FOC"], C["Z_FOC"], lam, C["ca_f"])
        _cache[key] = (P, D, ok & ok2)
    return _cache[key]

# ---------------- slit -> grating frame ----------------
a = np.radians(20.0); A = np.array([0.0, 96.9 * np.sin(a), 96.9 * np.cos(a)])
U_AX = -A / np.linalg.norm(A); X_H = np.array([1.0, 0, 0]); E_W = np.cross(X_H, U_AX)
Oc = np.array([0.0, 0.0, V.RSUB])

def at_plane(P, D, zs):
    return P + ((zs - P[:, 2]) / D[:, 2])[:, None] * D

def through_slit(P, D, ok, zs):
    Q = at_plane(P, D, zs)
    keep = ok & (np.abs(Q[:, 1]) <= SLIT_W / 2) & (np.abs(Q[:, 0]) <= SLIT_H / 2)
    return Q[keep], D[keep]

def to_grating(Q, D):
    S = A + np.outer(Q[:, 0], X_H) + np.outer(Q[:, 1], E_W)
    d = np.outer(D[:, 0], X_H) + np.outer(D[:, 1], E_W) + np.outer(D[:, 2], U_AX)
    return S, d / np.linalg.norm(d, axis=1, keepdims=True)

def grating_trace(S, d, lam, stop):
    e = S - Oc; ed = np.sum(d * e, 1); t = -ed + np.sqrt(np.clip(ed * ed - (np.sum(e * e, 1) - V.RSUB ** 2), 0, None))
    P = S + t[:, None] * d
    inside = (np.hypot(P[:, 0], P[:, 1]) <= EPD10 / 2) if stop else ((np.abs(P[:, 0]) <= 20) & (np.abs(P[:, 1]) <= 20))
    Nrm = (Oc - P) / V.RSUB; g = np.array([0, 1.0, 0]) - Nrm[:, 1:2] * Nrm; g /= np.linalg.norm(g, axis=1, keepdims=True)
    nu = R.nu_of(P[:, 1], 400.0, *CHIRP)
    din = d - np.sum(d * Nrm, 1, keepdims=True) * Nrm; dt = din + (lam * 1e-3 * nu)[:, None] * g
    n2 = 1 - np.sum(dt * dt, 1); ok = inside & (n2 >= 0)
    dout = dt + (-np.sign(np.sum(d * Nrm, 1)) * np.sqrt(np.clip(n2, 0, None)))[:, None] * Nrm
    Cd, ndet, ed_, ex = V.CS; H = P + (((Cd - P) @ ndet) / (dout @ ndet))[:, None] * dout
    return (H - Cd) @ ed_, (H - Cd) @ ex, ok

def plate_scale(l):
    b = np.arcsin(l * 1e-3 * 400 - np.sin(a)); r = V.plane_r(l)
    return (1e6 / 400) * np.cos(b) * np.cos(b - np.radians(V.PBETA)) / r, r   # Shimadzu eq (24): D = cos(b)cos(bi-b)/(N m r')

def fwhm(u_um):
    """FWHM of a traced line profile after convolving with the 8 um pixel (um)."""
    h, edges = np.histogram(u_um, bins=np.arange(-120, 120.25, 0.25)); h = np.convolve(h, np.ones(int(PIXW / 0.25)), "same")
    x = 0.5 * (edges[1:] + edges[:-1]); half = h.max() / 2; above = np.where(h >= half)[0]; i0, i1 = above[0], above[-1]
    xl = np.interp(half, [h[i0 - 1], h[i0]], [x[i0 - 1], x[i0]]); xr = np.interp(half, [h[i1 + 1], h[i1]], [x[i1 + 1], x[i1]])
    return xr - xl

def evaluate(C, zs, stop=False, detail=False):
    out = []
    for l in BAND:
        P, D, ok = relay_out(C, l); Q, Dq = through_slit(P, D, ok, zs)
        S, d = to_grating(Q, Dq); U, Vv, gok = grating_trace(S, d, l, stop)
        frac = gok.sum() / NRAY; Vc = np.median(Vv[gok]); onrow = gok & (np.abs(Vv - Vc) <= PIXH / 2e3)
        row = dict(lam=int(l * 1e3), relay_slit=float(frac), on_pixel=float(onrow.sum() / max(gok.sum(), 1)),
                   system=float(onrow.sum() / NRAY))
        if detail:
            Uc, _, okc = grating_trace(np.repeat(A[None, :], len(d), 0), d, l, stop)
            rms = (Uc[okc] - Uc[okc].mean()).std() * 1e3
            sinmax = np.percentile(np.linalg.norm(np.cross(d[gok], U_AX), axis=1), 99)
            feff = 10.0 if stop else 1 / (2 * sinmax); plate, r = plate_scale(l); diff = 1.028 * l * feff
            row.update(rms_um=float(rms), f_eff=float(feff),
                       dlam=float(np.hypot(np.hypot(SLIT_W * 1e3, 2.355 * rms), np.hypot(PIXW, diff)) / 1e3 * plate),
                       dlam_fwhm=float(np.hypot(fwhm((U[onrow] - np.median(U[onrow])) * 1e3), diff) / 1e3 * plate),
                       floor=float(np.hypot(SLIT_W * 1e3, np.hypot(PIXW, diff)) / 1e3 * plate))
            row["_U"] = ((U[gok] - np.median(U[onrow])) * 1e3)[:25000]; row["_V"] = ((Vv[gok] - Vc) * 1e3)[:25000]
            row["_Urow"] = (U[onrow] - np.median(U[onrow])) * 1e3
            if abs(l - 0.7) < 1e-9: row["_slit"] = at_plane(P, D, zs)[ok][:40000]
        out.append(row)
    return out

def ideal_reference(C, stop=False):
    NAo = NAF / C["m"]; D = C["m"] * CORE * 1e3; rr = D / 2; w = SLIT_W * 1e3 / 2
    frac = 2 * (w * np.sqrt(rr * rr - w * w) + rr * rr * np.arcsin(w / rr)) / (np.pi * rr * rr)
    return frac * (min(1.0, (0.05 / NAo) ** 2) if stop else 1.0)

def run(name):
    C = config(name)
    print(f"\n==== {name}: {C['label']}  m = {C['m']:.2f} ====")
    print(f"layout: collimator first vertex {C['D1']:.2f} mm | gap {C['gap']} mm | focusing lens at {C['Z_FOC']:.2f} mm | nominal slit {C['Z_SLIT0']:.2f} mm")
    zbest = minimize_scalar(lambda zs: -np.mean([r["system"] for r in evaluate(C, zs)]),
                            bounds=(C["Z_SLIT0"] - 4, C["Z_SLIT0"] + 4), method="bounded", options=dict(xatol=0.01)).x
    print(f"best slit position: {zbest:.2f} mm from the fiber face ({zbest - C['Z_SLIT0']:+.2f} mm from nominal)")
    res = evaluate(C, zbest, stop=False, detail=True); stp = evaluate(C, zbest, stop=True, detail=True)
    print(f"{'lam':>5} {'relay+slit':>11} {'on-pixel':>9} {'system':>8} {'dlam':>6} {'fwhm':>6} {'f/#':>5} | f/10 stop: {'system':>7} {'dlam':>6}")
    for r, q in zip(res, stp):
        print(f"{r['lam']:5d} {r['relay_slit']*100:10.2f}% {r['on_pixel']*100:8.1f}% {r['system']*100:7.3f}% {r['dlam']:6.2f} "
              f"{r['dlam_fwhm']:6.2f} {r['f_eff']:5.1f} | {q['system']*100:7.3f}% {q['dlam']:6.2f}")
    mean = lambda rows, k: float(np.mean([r[k] for r in rows]))
    m0 = mean(res, "system")
    tol = [(dz, mean(evaluate(C, zbest + dz), "system") / m0) for dz in (-1.0, -0.5, 0.5, 1.0)]
    print(f"mean: system {m0*100:.3f}% (no stop) / {mean(stp,'system')*100:.3f}% (f/10 stop) | relay+slit {mean(res,'relay_slit')*100:.2f}% "
          f"vs ideal {ideal_reference(C)*100:.2f}% | dlam {mean(res,'dlam'):.2f} nm (formula) {mean(res,'dlam_fwhm'):.2f} nm (direct FWHM)")
    print("focus tolerance: " + "  ".join(f"{dz:+.1f}mm:{v:.2f}" for dz, v in tol))
    strip = lambda rows: [{k: v for k, v in r.items() if not k.startswith("_")} for r in rows]
    json.dump(dict(pair=C["label"], magnification=C["m"],
                   layout=dict(d1=C["D1"], gap=C["gap"], z_foc=C["Z_FOC"], z_slit_nominal=C["Z_SLIT0"], z_slit=zbest),
                   no_stop=strip(res), f10_stop=strip(stp), ideal_relay_slit=ideal_reference(C), focus_tolerance=tol),
              open(f"relay_design_{name}.json", "w"), indent=1)
    performance_figure(C, res, stp, zbest); spots_figure(C, res)
    _cache.clear()
    return C, zbest

# ---------------- figures ----------------
def draw_layout(ax, C, zbest, y0=0.0):
    N = 9
    for l, col in zip(BAND, COLS):
        k = np.linspace(-NAF, NAF, N); P = np.zeros((N, 3)); D = np.stack([np.zeros(N), k, np.sqrt(1 - k * k)], 1)
        segs = [P.copy()]; Pp, Dp = P, D
        for surfs, z0, ca in ((C["COLL"], C["D1"], C["ca_c"]), (C["FOC"], C["Z_FOC"], C["ca_f"])):
            zz = z0
            for s in surfs:
                Pp, Dp, _, _ = RO.trace(Pp, Dp, [s], zz, l, ca); segs.append(Pp.copy()); zz += s[1]
        segs.append(at_plane(Pp, Dp, zbest))
        for i in range(N): ax.plot([sg[i, 2] for sg in segs], [y0 + sg[i, 1] for sg in segs], color=col, lw=0.6, alpha=0.8)
    for surfs, z0, rad in ((C["COLL"], C["D1"], C["rad_c"]), (C["FOC"], C["Z_FOC"], C["rad_f"])):
        zz = z0; ends = []
        for Rr, t, g in surfs:
            yy = np.linspace(-rad, rad, 100); sag = Rr - np.sign(Rr) * np.sqrt(Rr * Rr - yy * yy)
            ax.plot(zz + sag, y0 + yy, "k", lw=1); ends.append(zz + sag[-1]); zz += t
        for e0, e1 in zip(ends[:-1], ends[1:]):
            for sgn in (-1, 1): ax.plot([e0, e1], [y0 + sgn * rad, y0 + sgn * rad], "k", lw=1)
    ax.plot([0, 0], [y0 - 0.1, y0 + 0.1], color="#2f8f5b", lw=4)
    ax.plot([zbest, zbest], [y0 - 1.5, y0 + 1.5], color="k", lw=2.5)
    ax.annotate("", xy=(zbest, y0 - C["rad_f"] - 3), xytext=(0, y0 - C["rad_f"] - 3),
                arrowprops=dict(arrowstyle="<->", color="0.4", lw=0.8))
    ax.text(zbest / 2, y0 - C["rad_f"] - 6.5, f"fiber → slit {zbest:.0f} mm", ha="center", fontsize=9, color="0.25")
    ax.text(-4, y0 + C["rad_c"] + 2, C["label"], fontsize=9, fontweight="bold", va="bottom")

def layouts_figure(items):
    fig, ax = plt.subplots(figsize=(13, 6.2)); y = 0.0; ys = []
    for C, zb in items:
        draw_layout(ax, C, zb, y); ys.append(y); y -= 36
    zmax = max(zb for _, zb in items)
    ax.plot([0, zmax], [0, 0], alpha=0)
    ax.set_aspect("equal"); ax.set_xlim(-6, zmax + 8); ax.set_ylim(y + 17, 20); ax.set_yticks([]); ax.set_xlabel("mm (to scale)")
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    yb = y + 22; ax.plot([0, 96.9], [yb, yb], color="#1aa5cf", lw=3, alpha=0.7)
    ax.text(98.5, yb, "for scale: P0400 slit → grating arm, 96.9 mm", va="center", fontsize=9, color="#0c5468")
    ax.set_title("Fiber-to-slit relay, both catalog options to the same scale — real surfaces, 400–1100 nm rays", fontsize=11)
    fig.tight_layout(); fig.savefig("relay_layouts.png", dpi=160, bbox_inches="tight"); plt.close(fig)

def performance_figure(C, res, stp, zbest):
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(16, 4.5))
    Qs = [r for r in res if r["lam"] == 700][0]["_slit"]
    a1.scatter(Qs[:, 1] * 1e3, Qs[:, 0] * 1e3, s=0.3, color="#e8912a", alpha=0.35)
    a1.add_patch(plt.Rectangle((-SLIT_W * 500, -1500), SLIT_W * 1e3, 3000, fc="none", ec="k", lw=1.4))
    ri = C["m"] * CORE * 500; a1.add_patch(plt.Circle((0, 0), ri, fc="none", ec="#0c5468", ls="--", lw=1.2))
    a1.text(0, ri + 60, f"ideal image Ø{2*ri:.0f} µm", ha="center", fontsize=8, color="#0c5468"); a1.text(20, -1150, "25 µm slit", fontsize=8)
    a1.set_xlim(-1300, 1300); a1.set_ylim(-1300, 1300); a1.set_aspect("equal")
    a1.set_xlabel("across slit (µm)"); a1.set_ylabel("along slit (µm)"); a1.set_title("Fiber image at the slit, 700 nm", fontsize=10.5)
    a2.plot(LAM, [r["relay_slit"] * 100 for r in res], "o-", color="#1aa5cf", label="relay + slit")
    a2.plot(LAM, [r["system"] * 100 for r in res], "s-", color="#2f8f5b", label="system, no stop (recommended)")
    a2.plot(LAM, [r["system"] * 100 for r in stp], "s--", color="#e8912a", label="system, f/10 grating stop")
    a2.axhline(ideal_reference(C) * 100, color="#1aa5cf", ls=":", lw=1.2, label="ideal NA-matched relay")
    a2.set_xlabel("wavelength (nm)"); a2.set_ylabel("% of the fiber's light"); a2.set_title("Geometric throughput", fontsize=10.5)
    a2.grid(alpha=0.3); a2.legend(fontsize=8); a2.set_ylim(0, None)
    a3.plot(LAM, [r["dlam"] for r in res], "s-", color="#2f8f5b", label="no stop")
    a3.plot(LAM, [r["dlam"] for r in stp], "o--", color="#e8912a", label="f/10 grating stop")
    a3.set_xlabel("wavelength (nm)"); a3.set_ylabel("resolution Δλ (nm)"); a3.set_ylim(0, 2.5); a3.grid(alpha=0.3); a3.legend(fontsize=8)
    a3.set_title(f"End-to-end resolution — mean {np.mean([r['dlam'] for r in res]):.2f} nm", fontsize=10.5)
    fig.suptitle(f"Fiber #23-692 → {C['label']} → 25 µm slit #58-541 → P0400-01 → 8 × 200 µm pixels", fontsize=12, fontweight="bold", y=1.02)
    fig.tight_layout(); fig.savefig(f"relay_performance_{C['name']}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

def spots_figure(C, res):
    fig = plt.figure(figsize=(17, 10.2))
    gs = fig.add_gridspec(3, len(BAND), height_ratios=[1.35, 0.55, 1.0], hspace=0.28, wspace=0.12)
    bins = np.arange(-90, 90.5, 0.5); xc = 0.5 * (bins[1:] + bins[:-1])
    for i, (r, col) in enumerate(zip(res, COLS)):
        ax = fig.add_subplot(gs[0, i])
        ax.scatter(r["_U"], r["_V"], s=0.25, color=col, alpha=0.35, rasterized=True)
        ax.add_patch(plt.Rectangle((-PIXW / 2, -PIXH / 2), PIXW, PIXH, fc="none", ec="k", lw=1.0))
        ax.set_xlim(-80, 80); ax.set_ylim(-750, 750); ax.tick_params(labelsize=7.5); ax.set_xticklabels([])
        ax.set_title(f"{r['lam']} nm\nFWHM Δλ = {r['dlam_fwhm']:.2f} nm", fontsize=10)
        if i == 0: ax.set_ylabel("cross-dispersion (µm)", fontsize=9)
        else: ax.set_yticklabels([])
        ax.text(0.5, 0.02, f"{r['on_pixel']*100:.0f}% on pixel", transform=ax.transAxes, ha="center", fontsize=7.5, color="0.3")
        # line profile on the pixel row, convolved with the 8 um pixel
        bx = fig.add_subplot(gs[1, i]); h, _ = np.histogram(r["_Urow"], bins=bins); h = np.convolve(h, np.ones(int(PIXW / 0.5)), "same")
        h = h / h.max(); bx.fill_between(xc, h, color=col, alpha=0.55, lw=0); bx.plot(xc, h, color="k", lw=0.6)
        half = np.where(h >= 0.5)[0]; bx.plot([xc[half[0]], xc[half[-1]]], [0.5, 0.5], color="k", lw=1.2)
        bx.set_xlim(-80, 80); bx.set_ylim(0, 1.1); bx.set_yticks([]); bx.tick_params(labelsize=7.5)
        bx.set_xlabel("dispersion (µm)", fontsize=8)
        if i == 0: bx.set_ylabel("line\nprofile", fontsize=8.5)
    ax = fig.add_subplot(gs[2, 1:-1])
    ax.plot(LAM, [r["dlam_fwhm"] for r in res], "o-", color="#0c5468", lw=2.2,
            label=f"measured: FWHM of the traced line — mean {np.mean([r['dlam_fwhm'] for r in res]):.2f} nm")
    ax.plot(LAM, [r["dlam"] for r in res], "s--", color="#2f8f5b", lw=1.4,
            label=f"quick formula (Gaussian approximation) — mean {np.mean([r['dlam'] for r in res]):.2f} nm")
    ax.plot(LAM, [r["floor"] for r in res], ":", color="0.4", lw=1.4, label="floor: 25 µm slit + 8 µm pixel + diffraction")
    for r, col in zip(res, COLS): ax.plot(r["lam"], r["dlam_fwhm"], "o", color=col, ms=10, mec="k", mew=0.6, zorder=5)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("spectral resolution Δλ (nm)"); ax.set_ylim(0, 3.4); ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper center", ncol=1, frameon=False)
    fig.suptitle(f"End-to-end spot diagrams, line profiles and resolution — fiber 200 µm / NA 0.22 → {C['label']} → 25 µm slit → P0400-01\n"
                 "Dots: rays from the fiber that passed the slit. Box: one 8 × 200 µm pixel. Line profile: light on that pixel row "
                 "vs dispersion, with its FWHM (bar) = the resolution you would measure on a lamp line.", fontsize=10.5, y=0.995)
    fig.savefig(f"relay_spots_{C['name']}.png", dpi=150, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    names = [n for n in sys.argv[1:] if n in PAIRS] or list(PAIRS)
    items = [run(n) for n in names]
    if len(items) == len(PAIRS): layouts_figure(items)
