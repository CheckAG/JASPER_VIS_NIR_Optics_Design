#!/usr/bin/env python3
"""JASPER x Shimadzu P0400-01: the complete optical layout, to scale, with traced rays.

fiber (Edmund #23-692) -> compact relay (#49-321 reversed + #49-329) -> 25 um slit (#58-541)
-> P0400-01 concave flat-field grating -> flat detector on the datasheet plane.

Drawn in the grating's dispersion plane (horizontal = along the grating normal, vertical = dispersion
direction). Every surface is real: lens prescriptions from relay_optics.py, grating sphere R = 99.61 mm
with the solved groove chirp, detector plane from the datasheet. Relay rays are traced surface by surface
(at 700 nm -- the relay's colour error is invisible at this scale); after the slit the white-light beam is
launched from the slit centre with the relay's output cone, then diffracted wavelength by wavelength.
Needs relay_design_compact.json (python relay_design.py compact) and validate_p0400.json.
"""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon
import relay_optics as RO, relay_design as RD, validate_p0400 as V, vls_raytrace as R

C = RD.config("compact")
LAYOUT = json.load(open("relay_design_compact.json")); ZS = LAYOUT["layout"]["z_slit"]
PERF = LAYOUT["no_stop"]
A, U, EW = RD.A, RD.U_AX, RD.E_W
RED = "#d81f1f"; GLASS = "#dcecf5"; BAND = V.BAND
COLS = plt.cm.turbo(np.linspace(0.08, 0.92, len(BAND)))

def g(zr, yr):
    """relay coordinates (axial zr from the fiber face, transverse yr) -> plot coordinates (z, y)."""
    zr = np.atleast_1d(np.asarray(zr, float)); yr = np.atleast_1d(np.asarray(yr, float))
    P = A[None, :] + (zr - ZS)[:, None] * U[None, :] + yr[:, None] * EW[None, :]
    return P[:, 2], P[:, 1]

def surface_curve(Rr, zv, rad, n=60):
    yy = np.linspace(-rad, rad, n); return zv + (Rr - np.sign(Rr) * np.sqrt(Rr * Rr - yy * yy)), yy

def draw_lens(ax, surfs, z0, rad):
    zz = z0; curves = []
    for Rr, t, gl in surfs:
        curves.append(surface_curve(Rr, zz, rad)); zz += t
    for (z1, y1), (z2, y2) in zip(curves[:-1], curves[1:]):          # one filled polygon per element
        pz, py = g(np.concatenate([z1, z2[::-1]]), np.concatenate([y1, y2[::-1]]))
        ax.add_patch(Polygon(np.column_stack([pz, py]), closed=True, fc=GLASS, ec="k", lw=0.9, zorder=4))

def relay_fan(n=11, lam=0.7):
    k = np.linspace(-RD.NAF, RD.NAF, n); P = np.zeros((n, 3)); D = np.stack([np.zeros(n), k, np.sqrt(1 - k * k)], 1)
    pts = [P.copy()]
    for surfs, z0, ca in ((C["COLL"], C["D1"], C["ca_c"]), (C["FOC"], C["Z_FOC"], C["ca_f"])):
        zz = z0
        for s in surfs:
            P, D, _, _ = RO.trace(P, D, [s], zz, lam, ca); pts.append(P.copy()); zz += s[1]
    pts.append(RD.at_plane(P, D, ZS)); return pts, D

def grating_rays(dirs, lam):
    """launch from the slit centre with the given grating-frame directions; return grating hits P and detector hits H."""
    S = np.repeat(A[None, :], len(dirs), 0); Oc = RD.Oc; Rs = V.RSUB
    e = S - Oc; ed = np.sum(dirs * e, 1); t = -ed + np.sqrt(np.clip(ed * ed - (np.sum(e * e, 1) - Rs * Rs), 0, None))
    P = S + t[:, None] * dirs; Nrm = (Oc - P) / Rs
    gg = np.array([0, 1.0, 0]) - Nrm[:, 1:2] * Nrm; gg /= np.linalg.norm(gg, axis=1, keepdims=True)
    nu = R.nu_of(P[:, 1], 400.0, *RD.CHIRP)
    din = dirs - np.sum(dirs * Nrm, 1, keepdims=True) * Nrm; dt = din + (lam * 1e-3 * nu)[:, None] * gg
    dout = dt + (-np.sign(np.sum(dirs * Nrm, 1)) * np.sqrt(np.clip(1 - np.sum(dt * dt, 1), 0, None)))[:, None] * Nrm
    Cd, ndet, _, _ = V.CS; H = P + (((Cd - P) @ ndet) / (dout @ ndet))[:, None] * dout
    return P, H

def main():
    fig, ax = plt.subplots(figsize=(15, 8.6))
    # ---------------- grating blank (cross-section through the meridian) ----------------
    Rs = V.RSUB; yy = np.linspace(-20, 20, 200); zsurf = Rs - np.sqrt(Rs * Rs - yy * yy)
    ax.fill(np.concatenate([zsurf, [-9.5, -9.5]]), np.concatenate([yy, [20, -20]]), color="#cfe6f2", ec="k", lw=1, zorder=2)
    ax.plot(zsurf, yy, "k", lw=2, zorder=3)
    for yt in np.linspace(-19, 19, 20):
        zt = Rs - np.sqrt(Rs * Rs - yt * yt); ax.plot([zt - 0.5, zt + 0.5], [yt, yt], color="#3a6ea5", lw=0.5, zorder=3)
    ax.plot([0, 60], [0, 0], ls=(0, (7, 5)), color="0.55", lw=1); ax.text(61, 0, "grating normal", fontsize=8.5, va="center", color="0.4")
    ax.text(-11, -21.5, "Shimadzu P0400-01\n400 l/mm at the pole, R = 99.6 mm\n40 × 40 mm blank (spherical)", fontsize=8.5, ha="right", va="top")
    # ---------------- relay: rays, lenses, fiber, slit ----------------
    pts, Dout = relay_fan()
    for i in range(len(pts[0])):
        zc, yc = zip(*[g(p[i, 2], p[i, 1]) for p in pts])
        ax.plot(np.ravel(zc), np.ravel(yc), color=RED, lw=0.55, alpha=0.75, zorder=1)
    draw_lens(ax, C["COLL"], C["D1"], C["rad_c"]); draw_lens(ax, C["FOC"], C["Z_FOC"], C["rad_f"])
    fz, fy = g([-14, 0, 0, -14, -14], [-1.6, -1.6, 1.6, 1.6, -1.6])                     # SMA905 ferrule (3.18 mm)
    ax.add_patch(Polygon(np.column_stack([fz, fy]), closed=True, fc="0.75", ec="k", lw=0.8, zorder=4))
    for yr0, yr1 in ((0.0125, 2.2), (-0.0125, -2.2)):                                   # slit jaws
        jz, jy = g([ZS, ZS], [yr0, yr1]); ax.plot(jz, jy, color="k", lw=3.2, solid_capstyle="butt", zorder=6)
    mz, my = g([ZS - 1.25, ZS, ZS, ZS - 1.25, ZS - 1.25], [-12.5, -12.5, 12.5, 12.5, -12.5])   # Ø25 mm slit mount
    ax.add_patch(Polygon(np.column_stack([mz, my]), closed=True, fc="0.88", ec="0.45", lw=0.7, zorder=3))
    # ---------------- entrance beam slit -> grating, then diffracted rays ----------------
    dirs = np.outer(Dout[:, 1], EW) + np.outer(Dout[:, 2], U); dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    foci = []
    for l, col in zip(BAND, COLS):
        P, H = grating_rays(dirs, l)
        if l == BAND[0]:
            for p in P: ax.plot([A[2], p[2]], [A[1], p[1]], color=RED, lw=0.55, alpha=0.75, zorder=1)
        for p, h in zip(P, H): ax.plot([p[2], h[2]], [p[1], h[1]], color=col, lw=0.5, alpha=0.6, zorder=1)
        foci.append(H.mean(0))
    # detector: datasheet plane between the 400 and 1100 nm points (plus a little margin)
    e1, e2 = foci[0], foci[-1]; ev = (e2 - e1) / np.linalg.norm(e2 - e1)
    d1, d2 = e1 - 1.0 * ev, e2 + 1.0 * ev
    ax.plot([d1[2], d2[2]], [d1[1], d2[1]], color="k", lw=4, solid_capstyle="butt", zorder=6)
    for f, col, l in zip(foci, COLS, BAND):
        ax.plot(f[2], f[1], "o", color=col, ms=5.5, mec="k", mew=0.4, zorder=7)
        ax.text(f[2] + 2.2, f[1], f"{int(l*1e3)}", fontsize=7.5, va="center", color="0.25", zorder=7)
    ax.text(e1[2] - 2, e1[1] - 3.0, f"flat detector, {np.linalg.norm(e2-e1):.2f} mm (datasheet plane)", fontsize=8.5, va="top", ha="right")
    # ---------------- dimensions / labels ----------------
    mid = 0.5 * A; ax.text(mid[2] - 6, mid[1] + 3, "slit → grating 96.9 mm", fontsize=8.5, rotation=20, color=RED, style="italic")
    ax.add_patch(Arc((0, 0), 34, 34, theta1=0, theta2=20, lw=1)); ax.text(18.5, 2.3, "α = 20°", fontsize=9)
    fz0, fy0 = g(0, 0); lz, ly = g(ZS / 2, 9.5)
    ax.text(lz[0], ly[0], f"fiber → slit {ZS:.0f} mm", fontsize=8.5, rotation=20, ha="center", color="0.25")
    tz, ty = g(-14, 3.0); ax.text(tz[0] - 2, ty[0] + 1, "fiber #23-692\n200 µm, NA 0.22, SMA905", fontsize=8.5, ha="right")
    cz, cy = g(C["D1"] + 3, C["rad_c"] + 1.5); ax.text(cz[0], cy[0], "#49-321\n14 mm (reversed)", fontsize=8, ha="center")
    kz, ky = g(C["Z_FOC"] + 2.5, C["rad_f"] + 1.5); ax.text(kz[0], ky[0], "#49-329\n60 mm", fontsize=8, ha="center")
    sz, sy = g(ZS, -13.5); ax.text(sz[0] + 3, sy[0] + 1, "25 µm slit #58-541\n(Ø25 mm mount)", fontsize=8, ha="left", va="top")
    # clearance between slit mount and detector active end (a real mechanical check)
    mount_low = np.array(g(ZS, -12.5)).ravel(); top = e2 if e2[1] > e1[1] else e1
    clr = mount_low[1] - (top[1] + 1.0)
    # ---------------- info box ----------------
    mean_fwhm = np.mean([r["dlam_fwhm"] for r in PERF]); best = min(r["dlam_fwhm"] for r in PERF)
    light = np.mean([r["system"] for r in PERF]) * 100
    box = ("Fiber      #23-692  200 µm core, NA 0.22, SMA905\n"
           f"Relay      #49-321 14 mm (rev.) + #49-329 60 mm   ×{C['m']:.2f} → f/{1/(2*RD.NAF/C['m']):.1f}\n"
           "Slit       #58-541  25 µm × 3 mm\n"
           "Grating    P0400-01  400 l/mm, R 99.6 mm, α 20°\n"
           "Detector   flat, 400–1100 nm, 8 × 200 µm pixels\n"
           f"Resolution {mean_fwhm:.2f} nm mean FWHM ({best:.2f} nm at 500/1000 nm)\n"
           f"Light      ~{light:.2f} % of fiber light on the pixels\n"
           f"Clearance  slit mount ↔ detector active end ≈ {clr:.0f} mm")
    ax.text(0.985, 0.03, box, transform=ax.transAxes, fontsize=8.3, family="monospace", ha="right", va="bottom", ma="left",
            bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="0.6", lw=0.8))
    ax.set_aspect("equal"); ax.set_xlim(-45, 190); ax.set_ylim(-40, 78)
    ax.set_xlabel("mm (to scale) — along the grating normal"); ax.set_ylabel("mm — dispersion direction")
    ax.set_title("JASPER × Shimadzu P0400-01 — fiber, relay and flat-field grating, to scale, traced rays", fontsize=13)
    fig.text(0.5, 0.005, "All surfaces real (Edmund lens prescriptions, grating sphere with the solved groove chirp, datasheet detector plane). "
             "Relay rays at 700 nm; after the slit the beam is drawn from the slit centre with the relay's output cone. Modelled, not built.",
             ha="center", fontsize=8, color="0.35")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig("jasper_p0400_layout.png", dpi=160, bbox_inches="tight", facecolor="white")
    print(f"saved jasper_p0400_layout.png | slit-mount to detector clearance {clr:.1f} mm | detector {np.linalg.norm(e2-e1):.2f} mm")

if __name__ == "__main__":
    main()
