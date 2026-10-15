#!/usr/bin/env python3
"""Resolution metrics re-check: quick formula vs measured line width, for every headline design.

Why this exists. Earlier logs quoted resolution from a quick formula,
    dlam = hypot(slit, 2.355*RMS, pixel, diffraction) x plate scale,
which converts the RMS ray spread to a FWHM as if every line were bell-shaped (Gaussian). Real lines
are not, and the error goes BOTH ways depending on which aberration dominates:
  - defocus-dominated lines are FLAT-TOPPED -> the formula reads too narrow (optimistic)
  - coma-dominated lines have a sharp core and a long one-sided TAIL -> the formula reads wider than the
    core, but the FWHM alone ignores the tail, which can hold half the light
So we report two measured widths from the same traced line (finite slit traced as 41 points across its
width, 8 um pixel applied, diffraction added in quadrature, x plate scale from traced centroids):
  - FWHM: what you read off a lamp line
  - EE80: the narrowest spectral window holding 80 % of the line's energy -- catches tails the FWHM misses
and the old formula alongside. Where FWHM and EE80 disagree, the line has tails or a flat top.

Cases: Release 1 panels A/B/C (raytrace_three_panels.py), the realizable 488 nm two-point recording
(recording_optimize_sim.py), the same VLS toroid at f/5.3, and the Release 3 P0400 build (relay_design.py;
needs validate_p0400.json from step R1 in the working directory).
Writes jasper_resolution_metrics.png and resolution_metrics.json. Runtime ~6-8 min (two VLS fits).
"""
import numpy as np, json, os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vls_raytrace as R

SLIT = 20.0; PIX = 8.0
BAND = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]; LAM_NM = np.array(BAND) * 1e3
N = 230.0; ALPHA = 12.34; rA = 137.357; rB = 126.75; Ry = 133.0; a = np.radians(ALPHA)
_bc = np.arcsin(0.75e3 / (1e6 / N) - np.sin(a)); RHO = (np.cos(a) + np.cos(_bc)) / (1 / rA + 1 / rB)
PERP = np.array([0.0, np.cos(a), -np.sin(a)])          # slit-width direction (perpendicular to the incident chief ray)
SLIT_PTS = np.linspace(-SLIT / 2e3, SLIT / 2e3, 41)    # mm

# ---------------- line-width measures (u in um on the detector) ----------------
def pixel_profile(u, step=0.1):
    h, e = np.histogram(u, bins=np.arange(u.min() - 30, u.max() + 30, step)); k = int(round(PIX / step))
    return 0.5 * (e[1:] + e[:-1]), np.convolve(h, np.ones(k) / k, "same")

def fwhm(u):
    x, h = pixel_profile(u); half = h.max() / 2; ab = np.where(h >= half)[0]; i0, i1 = ab[0], ab[-1]
    return (np.interp(half, [h[i1 + 1], h[i1]], [x[i1 + 1], x[i1]]) - np.interp(half, [h[i0 - 1], h[i0]], [x[i0 - 1], x[i0]]))

def ee_width(u, frac=0.8):
    u = np.sort(u + np.random.default_rng(1).uniform(-PIX / 2, PIX / 2, len(u)))   # pixel as a uniform blur
    n = len(u); k = int(np.ceil(frac * n)); return float((u[k - 1:] - u[:n - k + 1]).min())

def measure(U_ext, rms_point, diff, plate):
    U = U_ext - np.median(U_ext)
    return dict(formula=float(np.hypot(np.hypot(SLIT, 2.355 * rms_point), np.hypot(PIX, diff)) / 1e3 * plate),
                fwhm=float(np.hypot(fwhm(U), diff) / 1e3 * plate),
                ee80=float(np.hypot(ee_width(U), diff) / 1e3 * plate)), U

# ---------------- simulated designs (vls_raytrace) ----------------
def sim_case(fn, rho, corr):
    EPD = rA / fn; d = R.Design("sim", N, ALPHA, rA, rB, Ry, rho, EPD, BAND)
    c = list(R.fit_vls(d)) if corr else [0.0, 0.0, 0.0]
    Cs = d.detector(c, d.best_focus(c)); A0 = d.A.copy()
    cen = [np.mean(d.trace(l, c, npup=1500, Csys=Cs)[0]) for l in BAND]; rows = []; prof = {}
    for i, l in enumerate(BAND):
        i0, i1 = max(0, i - 1), min(5, i + 1); plate = abs((LAM_NM[i1] - LAM_NM[i0]) / (cen[i1] - cen[i0]))
        U0, _ = d.trace(l, c, npup=2500, Csys=Cs); rms = (U0 - U0.mean()).std() * 1e3
        Us = []
        for s in SLIT_PTS:
            d.A = A0 + s * PERP; Us.append(d.trace(l, c, npup=1500, Csys=Cs)[0])
        d.A = A0.copy()
        m, U = measure(np.concatenate(Us) * 1e3, rms, 1.028 * l * (rB / EPD), plate)
        rows.append(dict(lam=int(l * 1e3), **m)); prof[int(l * 1e3)] = (U, plate)
    return dict(c=[float(x) for x in c]), rows, prof

def holo_case():
    import recording_optimize_sim as M
    j = json.load(open("recording_optimize_sim.json" if os.path.exists("recording_optimize_sim.json")
                       else os.path.join("figures", "recording_optimize_sim.json")))["best_recording"]
    C = M.src(j["rC"], j["gamma"]); D = M.src(j["rD"], j["delta"]); lam0 = j["lam0"] * 1e-6
    nu = lambda P, g: M.nu_holo(P, g, C, D, lam0)
    Cs = M.detector(j["rB"]); A0 = M.A.copy(); cen = [np.mean(M.trace(l, nu, Cs, 1200)[0]) for l in BAND]; rows = []
    for i, l in enumerate(BAND):
        i0, i1 = max(0, i - 1), min(5, i + 1); plate = abs((LAM_NM[i1] - LAM_NM[i0]) / (cen[i1] - cen[i0]))
        U0, _ = M.trace(l, nu, Cs, 2500); rms = (U0 - U0.mean()).std() * 1e3
        Us = []
        for s in SLIT_PTS:
            M.A = A0 + s * PERP; Us.append(M.trace(l, nu, Cs, 1500)[0])
        M.A = A0.copy()
        m, _ = measure(np.concatenate(Us) * 1e3, rms, 1.028 * l * (j["rB"] / M.EPD), plate); rows.append(dict(lam=int(l * 1e3), **m))
    return dict(lam0=j["lam0"]), rows, {}

def p0400_case():
    import relay_design as RD
    C = RD.config("compact")
    from scipy.optimize import minimize_scalar
    zs = minimize_scalar(lambda z: -np.mean([r["system"] for r in RD.evaluate(C, z)]),
                         bounds=(C["Z_SLIT0"] - 4, C["Z_SLIT0"] + 4), method="bounded", options=dict(xatol=0.01)).x
    rows = []; prof = {}
    for r in RD.evaluate(C, zs, detail=True):
        plate, _ = RD.plate_scale(r["lam"] / 1e3); diff = 1.028 * r["lam"] / 1e3 * r["f_eff"]; U = r["_Urow"]
        rows.append(dict(lam=r["lam"], formula=r["dlam"], fwhm=float(np.hypot(fwhm(U), diff) / 1e3 * plate),
                         ee80=float(np.hypot(ee_width(U), diff) / 1e3 * plate)))
        prof[r["lam"]] = (U - np.median(U), plate)
    return dict(slit_mm=float(zs)), rows, prof

CASES = [
    ("A", "Release 1 panel A\nf/5.3 spherical, classical", lambda: sim_case(5.3, Ry, False)),
    ("B", "Release 1 panel B\nf/8 toroidal, classical", lambda: sim_case(8.0, RHO, False)),
    ("C", "Release 1 panel C\nf/10 toroidal, VLS", lambda: sim_case(10.0, RHO, True)),
    ("H", "488 nm two-point\nrecording, f/10", holo_case),
    ("C53", "same VLS toroid\nat f/5.3", lambda: sim_case(5.3, RHO, True)),
    ("P", "P0400 build (R3)\nrelay + 25 µm slit", p0400_case),
]

def main():
    out = {}; profs = {}
    for tag, label, fn in CASES:
        meta, rows, prof = fn(); profs[tag] = prof
        mean = {k: float(np.mean([r[k] for r in rows])) for k in ("formula", "fwhm", "ee80")}
        worst = {k: float(max(r[k] for r in rows)) for k in ("formula", "fwhm", "ee80")}
        out[tag] = dict(label=label.replace("\n", " "), **meta, rows=rows, mean=mean, worst=worst)
        print(f"\n{label.replace(chr(10), ' | ')}")
        print("   lam    " + " ".join(f"{r['lam']:>6d}" for r in rows) + "    mean  worst")
        for k, name in (("formula", "formula"), ("fwhm", "FWHM"), ("ee80", "EE80")):
            print(f"   {name:8s}" + " ".join(f"{r[k]:6.2f}" for r in rows) + f"   {mean[k]:5.2f}  {worst[k]:5.2f}", flush=True)
    json.dump(out, open("resolution_metrics.json", "w"), indent=1)
    figure(out, profs)

COL = dict(formula="#2a78d6", fwhm="#eb6834", ee80="#1baf7a"); INK = "#0b0b0b"; INK2 = "#52514e"; GRID = "#e4e3df"

def figure(out, profs):
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig = plt.figure(figsize=(14, 9.2)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.42, wspace=0.22)
    shapes = [("A", 500, "Release 1 panel A (f/5.3, classical)", "Coma tail — FWHM sees only the core"),
              ("C", 700, "Release 1 panel C (f/10, VLS)", "Compact — all three measures agree"),
              ("P", 800, "P0400 build (Release 3)", "Flat top (defocus) — the formula reads narrow")]
    for k, (tag, lam, name, title) in enumerate(shapes):
        ax = fig.add_subplot(gs[0, k]); U, plate = profs[tag][lam]; x, h = pixel_profile(U); xn = x * plate / 1e3
        ax.fill_between(xn, h / h.max(), color="#c9c8c2", lw=0); ax.plot(xn, h / h.max(), color=INK2, lw=1)
        r = next(r for r in out[tag]["rows"] if r["lam"] == lam)
        # FWHM bar at half maximum, EE80 bracket below it
        half = np.where(h >= h.max() / 2)[0]; xl, xr = xn[half[0]], xn[half[-1]]
        ax.plot([xl, xr], [0.5, 0.5], color=COL["fwhm"], lw=2.5)
        us = np.sort(U); n = len(us); kk = int(np.ceil(0.8 * n)); j = int(np.argmin(us[kk - 1:] - us[:n - kk + 1]))
        el, er = us[j] * plate / 1e3, us[j + kk - 1] * plate / 1e3
        ax.plot([el, er], [0.12, 0.12], color=COL["ee80"], lw=2.5); ax.plot([el, el], [0.07, 0.17], color=COL["ee80"], lw=2.5)
        ax.plot([er, er], [0.07, 0.17], color=COL["ee80"], lw=2.5)
        ax.text(0.02, 0.97, f"formula  {r['formula']:.2f} nm\nFWHM     {r['fwhm']:.2f} nm\nEE80     {r['ee80']:.2f} nm",
                transform=ax.transAxes, va="top", ha="left", family="monospace", fontsize=8.5, color=INK)
        ax.set_title(f"{name} · {lam} nm\n{title}", fontsize=9.5, color=INK, loc="left")
        ax.set_xlim(-4.2, 2.2) if tag == "A" else ax.set_xlim(-2.2, 2.2); ax.set_ylim(0, 1.08)
        ax.set_xlabel("wavelength offset (nm)"); ax.set_yticks([]); ax.spines[["top", "right", "left"]].set_visible(False)
    ax = fig.add_subplot(gs[1, :]); tags = [t for t, _, _ in CASES]; xs = np.arange(len(tags)); w = 0.26
    for j, (k, name) in enumerate((("formula", "quick formula (2.355 × RMS, Gaussian)"), ("fwhm", "FWHM of the traced line"),
                                   ("ee80", "EE80: window holding 80 % of the light"))):
        v = [out[t]["mean"][k] for t in tags]; xb = xs + (j - 1) * w
        ax.bar(xb, v, width=w - 0.03, color=COL[k], label=name, zorder=3)
        for xi, vi in zip(xb, v): ax.text(xi, vi + 0.04, f"{vi:.2f}", ha="center", va="bottom", fontsize=8.5, color=INK)
    ax.set_xticks(xs); ax.set_xticklabels([lbl for _, lbl, _ in CASES], fontsize=8.8, color=INK)
    ax.set_ylabel("mean resolution over the band (nm)"); ax.set_ylim(0, max(out[t]["mean"]["formula"] for t in tags) * 1.18)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    ax.set_title("Mean resolution by three measures (500–1000 nm; P0400 400–1100 nm). Finite slit traced, 8 µm pixel, diffraction included.",
                 fontsize=9.5, color=INK, loc="left")
    fig.suptitle("JASPER resolution metrics re-check — the quick formula errs both ways; FWHM alone misses coma tails",
                 fontsize=12, color=INK, x=0.07, ha="left", y=0.995)
    fig.savefig("jasper_resolution_metrics.png", dpi=150, bbox_inches="tight"); print("\nsaved jasper_resolution_metrics.png and resolution_metrics.json")

if __name__ == "__main__":
    main()
