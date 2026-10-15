#!/usr/bin/env python3
"""Entrance options for the Shimadzu P0400-01 fed by an NA 0.22 fiber: resolution vs light on the pixels.

Grating model: validate_p0400.py (400 l/mm at the pole, alpha 20 deg, r 96.9 mm, R 99.61 mm from the blank,
detector on the datasheet plane, solved groove chirp per f-number, from validate_p0400.json).

Options
  A   fiber core IS the slit (fiber butted at the entrance, aperture mask at the grating sets the f-number)
  B   relay images the fiber onto a 20 um slit with magnification m = NA_fiber / NA_grating (NA-matched)
  B'  20 um slit butted directly to the fiber, no relay
Cores 50-400 um, f/5.3 and f/10.

Per option we compute
  resolution  -- traced aberration RMS per wavelength + slit (FWHM: 20 um, or 0.866*core for a round core)
                 + 8 um pixel + diffraction, in quadrature, times the plate scale (same formula as the rest
                 of the package)
  light       -- (fraction of the fiber's light accepted at the entrance: etendue)
                 x (fraction that lands on a 200 um-tall pixel row, traced from sources spread over the
                    entrance height -- P0400's spherical blank smears each wavelength into a tall astigmatic line)
                 x (fiber power ~ core area: an extended diffuse source fills the fiber)
                 normalised to option B, 100 um core, f/10.
Ideal relay here; relay_design.py traces the real one.

Run:  python p0400_entrance_options.py   -> p0400_entrance_options.json, p04_entrance_options.png
"""
import numpy as np, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import vls_raytrace as R
import validate_p0400 as V

NAF = 0.22; PIXH = 200.0; PIXW = 8.0; SLITW = 20.0
CHIRP = json.load(open("validate_p0400.json"))            # solved c1..c3 per f-number
CASES = [("A", 50), ("A", 100), ("B", 50), ("B", 100), ("B", 200), ("B", 400), ("B'", 100)]
FNUMS = [("f/5.3", 5.3), ("f/10", 10.0)]

def strip_frac(D, w):
    """fraction of a uniformly lit disc (diameter D) passed by a central slit of width w"""
    r = D / 2; a = min(w / 2, r)
    return 2 * (a * np.sqrt(r * r - a * a) + r * r * np.arcsin(a / r)) / (np.pi * r * r)

def vertical_capture(des, c, halfh, lam, ns=11, npup=500):
    """fraction landing within a 200 um-tall pixel row, for an entrance of half-height halfh (um), disc-weighted"""
    A0 = des.A.copy(); xs = np.linspace(-halfh, halfh, ns)
    wts = np.sqrt(np.clip(1 - (xs / halfh) ** 2, 0, None)) if halfh > 0 else np.ones(1)
    hits = []
    for x, wt in zip(xs, wts):
        des.A = A0 + np.array([x * 1e-3, 0, 0]); _, Vv = des.trace(lam, c, npup=npup, Csys=V.CS); hits.append((Vv, wt))
    des.A = A0
    Vc = np.median(np.concatenate([h[0] for h in hits]))
    num = sum(wt * np.sum(np.abs(Vv - Vc) * 1e3 <= PIXH / 2) for Vv, wt in hits)
    den = sum(wt * len(Vv) for Vv, wt in hits)
    return num / den

def main():
    rows = []
    for tag, fn in FNUMS:
        EPD = 96.9 / fn; NAS = 1 / (2 * fn); m = NAF / NAS; c = CHIRP[tag]["c"]
        des = R.Design("P0400", 400., 20., 96.9, 96.9, V.RSUB, V.RSUB, EPD, V.BAND)
        # aberration part of the resolution, once per f-number
        ab = []
        for l in V.BAND:
            U, _ = des.trace(l, c, npup=2500, Csys=V.CS); rms = (U - U.mean()).std() * 1e3
            b = np.arcsin(l * 1e-3 * 400 - np.sin(np.radians(20))); r = V.plane_r(l)
            plate = (1e6 / 400) * np.cos(b) * np.cos(b - np.radians(V.PBETA)) / r   # Shimadzu eq (24)
            ab.append((rms, plate, 1.028 * l * r / EPD))
        for opt, core in CASES:
            if opt == "A":  slit, acc, halfh = 0.866 * core, (NAS / NAF) ** 2, core / 2
            if opt == "B":  slit, acc, halfh = SLITW, strip_frac(m * core, SLITW), m * core / 2
            if opt == "B'": slit, acc, halfh = SLITW, strip_frac(core, SLITW) * (NAS / NAF) ** 2, core / 2
            dl = [np.hypot(np.hypot(slit, 2.355 * rms), np.hypot(PIXW, diff)) / 1e3 * plate for rms, plate, diff in ab]
            vcap = float(np.mean([vertical_capture(des, c, halfh, l) for l in (0.5, 0.75, 1.0)]))
            rows.append(dict(option=opt, core=core, fnum=tag, f=fn, slit_um=float(slit), acceptance=float(acc),
                             on_pixel=vcap, power=float(acc * vcap * (core / 100) ** 2),
                             mean=float(np.mean(dl)), best=float(min(dl)), worst=float(max(dl))))
    ref = [r for r in rows if r["option"] == "B" and r["core"] == 100 and r["fnum"] == "f/10"][0]["power"]
    for r in rows: r["light"] = r["power"] / ref
    json.dump(rows, open("p0400_entrance_options.json", "w"), indent=1)
    print(f"{'f#':6s} {'opt':3s} {'core':>5s} {'slit term':>10s} {'accept':>7s} {'on-pixel':>9s} {'light':>6s} {'mean':>6s} {'best':>5s} {'worst':>6s}")
    for r in rows:
        print(f"{r['fnum']:6s} {r['option']:3s} {r['core']:5d} {r['slit_um']*25.3/1e3:8.2f}nm {r['acceptance']*100:6.1f}% "
              f"{r['on_pixel']*100:8.1f}% {r['light']:5.2f}x {r['mean']:5.2f} {r['best']:5.2f} {r['worst']:6.2f}")
    make_figure(rows)

def make_figure(rows):
    COL = {5.3: "#e8912a", 10.0: "#1aa5cf"}; MK = {"A": "o", "B": "s", "B'": "^"}
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for r in rows:
        L, mn = r["light"], r["mean"]
        ax.errorbar(L, mn, yerr=[[mn - r["best"]], [r["worst"] - mn]], fmt=MK[r["option"]], color=COL[r["f"]], ms=9,
                    mec="k", mew=0.6, elinewidth=0.9, capsize=3)
        dx, dy, ha = 1.06, 0.07, "left"
        if r["option"] == "B'" and r["f"] == 5.3: dx, ha = 0.93, "right"
        if r["option"] == "B" and r["core"] == 50 and r["f"] == 5.3: dy = -0.22
        ax.text(L * dx, mn + dy, f"{r['core']} µm", fontsize=8, color="0.25", ha=ha)
    ax.set_xscale("log"); ax.set_xlim(0.12, 5.5); ax.set_ylim(0, 4.8)
    ax.set_xticks([0.125, 0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(["⅛×", "0.25×", "0.5×", "1×", "2×", "4×"]); ax.minorticks_off()
    ax.add_patch(FancyBboxPatch((0.88, 1.02), 2.15 - 0.88, 0.36, boxstyle="round,pad=0.02", fc="#e9f7fc", ec="#0c5468", lw=1.2, zorder=0))
    ax.text(1.38, 0.78, "recommended: B at f/10, 100–200 µm fiber", fontsize=9, color="#0c5468", fontweight="bold", ha="center")
    ax.set_xlabel("light reaching the pixels, relative to 100 µm fiber + relay + 20 µm slit at f/10")
    ax.set_ylabel("resolution Δλ, mean over 400–1100 nm (nm)\nbars: best–worst wavelength")
    ax.set_title("P0400-01 fed by an NA 0.22 fiber — entrance options (TCD1304-class 200 µm pixels)", fontsize=11.5, fontweight="bold")
    h = [Line2D([], [], marker="s", ls="", color="0.5", mec="k", label="B: relay (NA-matched) + 20 µm slit"),
         Line2D([], [], marker="o", ls="", color="0.5", mec="k", label="A: fiber core is the slit (no relay)"),
         Line2D([], [], marker="^", ls="", color="0.5", mec="k", label="B': 20 µm slit butted to fiber (no relay)"),
         Line2D([], [], marker="o", ls="", color=COL[5.3], label="f/5.3"), Line2D([], [], marker="o", ls="", color=COL[10.0], label="f/10")]
    ax.legend(handles=h, fontsize=8.5, loc="upper right", frameon=False); ax.grid(alpha=0.25)
    fig.text(0.5, -0.04, "Light assumes the fiber is filled by an extended diffuse source (e.g. a reflectance probe), so a bigger core carries more power. "
             "Relay treated as ideal here;\nits real colour error sets how evenly the slit is filled across 400–1100 nm (throughput), not the resolution.",
             ha="center", fontsize=8, color="0.35")
    fig.tight_layout(); fig.savefig("p04_entrance_options.png", dpi=150, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    main()
