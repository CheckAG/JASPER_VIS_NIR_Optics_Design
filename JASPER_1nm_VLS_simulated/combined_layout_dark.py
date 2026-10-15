#!/usr/bin/env python3
"""Dark (Beam4-style) render of the complete JASPER x P0400-01 layout, matching raytrace_layout.py.

Same traced geometry as combined_layout.py -- fiber #23-692 -> compact relay (#49-321 reversed +
#49-329) -> 25 um slit -> P0400-01 -> datasheet detector plane -- drawn with dense low-alpha fans so the
beams glow: white light through the relay and into the grating, wavelength-coloured after it.
Needs relay_design_compact.json (python relay_design.py compact) and validate_p0400.json.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Arc
import combined_layout as CL
import relay_design as RD, validate_p0400 as V

BG = "#0e1116"; GREEN = "#39d98a"; WHITE = "#eaf2f8"; LABEL = "#dfe7ee"; ARROW = "#9fb0bf"
BAND = V.BAND; COLS = plt.cm.turbo(np.linspace(0.12, 0.95, len(BAND)))
C, ZS, A = CL.C, CL.ZS, CL.A

def lens_dark(ax, surfs, z0, rad):
    zz = z0; curves = []
    for Rr, t, gl in surfs:
        curves.append(CL.surface_curve(Rr, zz, rad)); zz += t
    for (z1, y1), (z2, y2) in zip(curves[:-1], curves[1:]):
        pz, py = CL.g(np.concatenate([z1, z2[::-1]]), np.concatenate([y1, y2[::-1]]))
        ax.add_patch(Polygon(np.column_stack([pz, py]), closed=True, fc="#7fb3d5", alpha=0.16, ec="none", zorder=4))
        ax.add_patch(Polygon(np.column_stack([pz, py]), closed=True, fc="none", ec="#9cc9e6", lw=1.0, zorder=5))

def note(ax, text, xy, xytext, color=LABEL, ha="left"):
    ax.annotate(text, xy=xy, xytext=xytext, color=color, fontsize=9, ha=ha,
                arrowprops=dict(arrowstyle="->", color=ARROW if color == LABEL else color, lw=0.9))

def main():
    fig, ax = plt.subplots(figsize=(15, 8.8)); fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
    # grating surface (exact sphere, R = 99.61 mm) over the 40 mm blank, with groove ticks
    Rs = V.RSUB; yy = np.linspace(-20, 20, 240); zz = Rs - np.sqrt(Rs * Rs - yy * yy)
    ax.plot(zz, yy, color=GREEN, lw=3.0, zorder=6, solid_capstyle="round")
    for yt in np.linspace(-19, 19, 30):
        zt = Rs - np.sqrt(Rs * Rs - yt * yt); ax.plot([zt, zt + 1.1], [yt, yt], color=GREEN, lw=0.6, alpha=0.55, zorder=5)
    # relay: dense white fan from the fiber through the real lens surfaces to the slit
    pts, Dout = CL.relay_fan(n=160)
    for i in range(len(pts[0])):
        zc, yc = zip(*[CL.g(p[i, 2], p[i, 1]) for p in pts])
        ax.plot(np.ravel(zc), np.ravel(yc), color=WHITE, lw=0.35, alpha=0.09, zorder=2)
    lens_dark(ax, C["COLL"], C["D1"], C["rad_c"]); lens_dark(ax, C["FOC"], C["Z_FOC"], C["rad_f"])
    fz, fy = CL.g([-14, 0, 0, -14, -14], [-1.6, -1.6, 1.6, 1.6, -1.6])
    ax.add_patch(Polygon(np.column_stack([fz, fy]), closed=True, fc="#5a6672", ec="#9fb0bf", lw=0.8, zorder=5))
    mz, my = CL.g([ZS - 1.25, ZS, ZS, ZS - 1.25, ZS - 1.25], [-12.5, -12.5, 12.5, 12.5, -12.5])
    ax.add_patch(Polygon(np.column_stack([mz, my]), closed=True, fc="#2a323c", ec="#5a6672", lw=0.7, zorder=3))
    for yr0, yr1 in ((0.0125, 2.2), (-0.0125, -2.2)):
        jz, jy = CL.g([ZS, ZS], [yr0, yr1]); ax.plot(jz, jy, color="#c9d2da", lw=3.2, solid_capstyle="butt", zorder=7)
    # slit -> grating white beam, then wavelength fans to the detector
    dirs = np.outer(Dout[:, 1], CL.EW) + np.outer(Dout[:, 2], CL.U); dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    foci = []
    for k, (l, col) in enumerate(zip(BAND, COLS)):
        P, H = CL.grating_rays(dirs, l)
        if k == 0:
            for p in P: ax.plot([A[2], p[2]], [A[1], p[1]], color=WHITE, lw=0.35, alpha=0.07, zorder=2)
        for p, h in zip(P, H): ax.plot([p[2], h[2]], [p[1], h[1]], color=col, lw=0.35, alpha=0.18, zorder=3)
        foci.append(H.mean(0))
    e1, e2 = foci[0], foci[-1]; ev = (e2 - e1) / np.linalg.norm(e2 - e1); d1, d2 = e1 - ev, e2 + ev
    ax.plot([d1[2], d2[2]], [d1[1], d2[1]], color="#c9d2da", lw=5, solid_capstyle="butt", zorder=7)
    for f, col, l in zip(foci, COLS, BAND):
        ax.plot(f[2], f[1], "o", color=col, ms=5.5, mec="k", mew=0.4, zorder=9)
    ax.text(e1[2] - 2, e1[1] - 0.5, "400 nm", color=COLS[0], fontsize=8.5, ha="right", va="center")
    ax.text(e2[2] + 2, e2[1] + 0.5, "1100 nm", color=COLS[-1], fontsize=8.5, ha="left", va="center")
    # normal, alpha, labels
    ax.plot([0, 60], [0, 0], ls=(0, (6, 5)), color="#5a6672", lw=1.0, zorder=1)
    ax.text(61, 0, "grating normal", color="#9fb0bf", fontsize=8.5, va="center", ha="left", zorder=8,
            bbox=dict(boxstyle="round,pad=0.2", fc=BG, ec="none", alpha=0.85))
    ax.add_patch(Arc((0, 0), 34, 34, theta1=0, theta2=20, lw=0.9, color=ARROW)); ax.text(18.5, 2.3, "α = 20°", color=LABEL, fontsize=9)
    note(ax, "Shimadzu P0400-01 flat-field grating\n400 l/mm · R = 99.6 mm · spherical", (Rs - np.sqrt(Rs * Rs - 400), 20), (-38, 31), color=GREEN)
    fz0, fy0 = CL.g(-7, -1.6); note(ax, "Fiber #23-692\n200 µm · NA 0.22 · SMA905", (fz0[0], fy0[0]), (fz0[0] + 6, fy0[0] - 22), ha="right")
    cz, cy = CL.g(C["D1"] + 4, C["rad_c"]); note(ax, "#49-321  14 mm\n(reversed)", (cz[0], cy[0]), (cz[0] - 4, cy[0] + 10), ha="center")
    kz, ky = CL.g(C["Z_FOC"] + 2.5, C["rad_f"]); note(ax, "#49-329  60 mm", (kz[0], ky[0]), (kz[0] - 16, ky[0] + 9), ha="center")
    note(ax, "Entrance slit  25 µm × 3 mm\n(#58-541, Ø25 mm mount)", (A[2], A[1]), (A[2] + 26, A[1] - 16))
    mid = (e1 + e2) / 2; note(ax, f"Detector — flat, {np.linalg.norm(e2 - e1):.2f} mm (datasheet plane)", (mid[2], mid[1]), (mid[2] + 12, mid[1] - 11))
    lz, ly = CL.g(ZS * 0.62, 7.5); ax.text(lz[0], ly[0], f"relay  fiber → slit {ZS:.0f} mm", color="#7f8c99", fontsize=8.5, rotation=20, ha="center")
    ax.text(0.5 * A[2] - 4, 0.5 * A[1] + 3.5, "slit → grating 96.9 mm", color="#7f8c99", fontsize=8.5, rotation=20, style="italic")
    # info box
    PERF = CL.PERF; mean_fwhm = np.mean([r["dlam_fwhm"] for r in PERF]); best = min(r["dlam_fwhm"] for r in PERF)
    light = np.mean([r["system"] for r in PERF]) * 100
    box = (f"Relay      ×{C['m']:.2f}  →  f/{1/(2*RD.NAF/C['m']):.1f} at the grating\n"
           "Band       400–1100 nm on 8 × 200 µm pixels\n"
           f"Resolution {mean_fwhm:.2f} nm mean FWHM  ({best:.2f} nm at 500 / 1000 nm)\n"
           f"Light      ~{light:.2f} % of the fiber's light on the pixels")
    ax.text(0.985, 0.04, box, transform=ax.transAxes, fontsize=8.6, family="monospace", ha="right", va="bottom", ma="left",
            color=LABEL, bbox=dict(boxstyle="round,pad=0.55", fc="#161b22", ec="#3a4450", lw=0.8))
    nrays = len(pts[0])
    ax.set_title(f"JASPER × Shimadzu P0400-01 — complete ray-traced optical layout   ·   {nrays} rays per wavelength, 8 wavelengths",
                 color="#eef3f7", fontsize=12.5, pad=14)
    ax.text(0.5, -0.03, "Real geometric ray trace, to scale, in the dispersion plane: fiber → real lens surfaces (Edmund prescriptions) → 25 µm slit → "
            "exact concave surface with the solved groove chirp → flat detector on the datasheet plane.  Relay rays at 700 nm; "
            "after the slit the beam is launched from the slit centre with the relay's output cone.  Modelled, not built.",
            transform=ax.transAxes, color=ARROW, fontsize=8.4, ha="center", va="top", wrap=True)
    ax.set_aspect("equal"); ax.set_xlim(-45, 196); ax.set_ylim(-40, 80); ax.axis("off")
    fig.tight_layout(); fig.savefig("jasper_p0400_layout_dark.png", dpi=160, bbox_inches="tight", facecolor=BG)
    print("wrote jasper_p0400_layout_dark.png")

if __name__ == "__main__":
    main()
