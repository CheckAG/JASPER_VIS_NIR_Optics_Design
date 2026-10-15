#!/usr/bin/env python3
"""Cross-validation: the analytic light-path-function (LPF) model vs the standalone ray tracer.

The machine-checkable form of "modelled two independent ways". For each design BOTH models are run on
the IDENTICAL geometry -- same grating, same band, and the SAME FIXED FLAT DETECTOR PLANE (normal along
beta_ref at perpendicular distance rB; the distance along the lambda ray is rB / cos(beta - beta_ref)).
Nothing is refocused. Metrics, with a delta column each:

  - plate scale (nm/mm)          -- grating equation; must agree to a fraction of a percent
  - astigmatic line height (mm)  -- F02 geometry; must agree closely (near zero for the toroid)
  - CLASSICAL resolution (nm)    -- constant grooves, no unknowns: the clean cross-check. Any residual
                                    gap is the higher-order aberration the 3rd-order LPF truncates;
                                    it grows with aberration strength (fast f, dense grooves).
  - CORRECTED resolution, split into two like-for-like cases (both models computed for each):
        IDEALISED  -- defocus nulled at EVERY wavelength. LPF: F20 = 0 at each lambda (coma only).
                      Ray trace: a separate best chirp fitted for each wavelength. This is a BOUND:
                      no real grating reaches it on a flat detector, because a single groove chirp
                      (VLS polynomial or two-point holographic recording) adds a defocus correction
                      ~ lambda x const and so nulls defocus at only ~2 wavelengths.
        ONE CHIRP  -- the physically real case: one linear chirp c1, IDENTICAL in both models
                      (fitted once in the ray tracer, then fed to the LPF as F20 + lambda*c1).
                      The ray tracer's full c1..c3 fit (c2 adds coma correction) is printed as well.
    Sign convention: the ray tracer measures y toward the entrance-slit side, the LPF's w the other
    way, so c1(LPF) = -c1(trace); the script checks this numerically at mid-band.

Designs: the simulated JASPER grating (spherical and toroidal; 500-1000 nm) and the real Shimadzu
P0400-01 using its datasheet geometry (R from the blank sag, detector on the datasheet plane, 400-1100 nm).

Run:  python validate_models.py        # table + jasper_model_validation.png + model_validation.json
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import lpf_resolution as L
import vls_raytrace as R

SLIT, PIX = 20.0, 8.0     # um -- same as lpf_resolution / vls_raytrace

# simulated toroid sag radius that nulls astigmatism at mid-band (same formula as lpf_resolution)
_a = np.radians(12.34); _rA = 137.357; _rB = 126.75
_bc = np.arcsin((0.75 * 1000) / (1e6 / 230.0) - np.sin(_a))
RHO_SIM = (np.cos(_a) + np.cos(_bc)) / (1.0 / _rA + 1.0 / _rB)

# P0400-01 datasheet geometry: R from the blank (9.5 mm centre, 13.6 mm at the 28.28 mm corner);
# detector plane from the datasheet mounting parameters: r'i = 96.93 mm (printed 96.9) at beta_i = -7.51 deg
R_P0400 = ((20 * np.sqrt(2)) ** 2 + 4.1 ** 2) / (2 * 4.1)     # 99.61 mm

SIM_BAND = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
P04_BAND = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
DESIGNS = [
    dict(label="Simulated VLS - SPHERICAL", N=230, alpha=12.34, rA=137.357, rB=126.75,
         Ry=133.0, rho=133.0, EPD=26.0, beta_ref=-2.4, band=SIM_BAND),
    dict(label="Simulated VLS - TOROIDAL", N=230, alpha=12.34, rA=137.357, rB=126.75,
         Ry=133.0, rho=RHO_SIM, EPD=26.0, beta_ref=-2.4, band=SIM_BAND),
    dict(label="Shimadzu P0400-01 (real)", N=400, alpha=20.0, rA=96.9, rB=96.93,
         Ry=R_P0400, rho=R_P0400, EPD=18.3, beta_ref=-7.51, band=P04_BAND),
]


def plane(d):
    """The fixed flat detector both models share: normal along beta_ref, perpendicular distance rB."""
    b = np.radians(d["beta_ref"])
    n = np.array([0.0, np.sin(b), np.cos(b)])
    ndet = -n; ex = np.array([1.0, 0, 0]); ed = np.cross(ndet, ex); ed /= np.linalg.norm(ed)
    return (d["rB"] * n, ndet, ed, ex)


def beta(d, l):
    return np.arcsin(l * 1e-3 * d["N"] - np.sin(np.radians(d["alpha"])))


def lpf_one_chirp(d, c1_trace):
    """LPF resolution with ONE linear chirp: F20_eff = F20 + lambda*c1_trace (see sign convention).
    Coma (F30) is left uncorrected, i.e. c2 = c3 = 0 -- matches the ray tracer's c1-only case."""
    a = np.radians(d["alpha"]); rA, Rt, EPD = d["rA"], d["Ry"], d["EPD"]; bref = np.radians(d["beta_ref"])
    w = np.linspace(-EPD / 2, EPD / 2, 201); out = []
    for l in d["band"]:
        b = beta(d, l); r = d["rB"] / np.cos(b - bref)
        F20 = np.cos(a)**2/rA - np.cos(a)/Rt + np.cos(b)**2/r - np.cos(b)/Rt
        F30 = (np.sin(a)/rA)*(np.cos(a)**2/rA - np.cos(a)/Rt) + (np.sin(b)/r)*(np.cos(b)**2/r - np.cos(b)/Rt)
        F2 = F20 + l * 1e-3 * c1_trace
        dy = (r / np.cos(b)) * (F2 * w + 1.5 * F30 * w**2); dy -= dy.mean()
        rms = np.sqrt(np.mean(dy**2)) * 1e3
        plate = (1e6 / d["N"]) * np.cos(b) * np.cos(b - bref) / r   # Shimadzu eq (24)
        out.append(np.hypot(np.hypot(SLIT, 2.355 * rms), np.hypot(PIX, 1.028 * l * r / EPD)) / 1e3 * plate)
    return float(np.mean(out))


def run_lpf(d):
    o = L.analyse(d["label"], d["N"], d["alpha"], d["rA"], d["rB"], d["Ry"], d["rho"],
                  d["EPD"], d["band"], d["beta_ref"])
    return dict(plate=float(np.mean([r["plate"] for r in o["rows"]])), astig=o["astig_mm"],
                classical=o["dlam_classical_mean"], idealised=o["dlam_corr_mean"])


def trace_eval(d, des, c, Cs):
    band = d["band"]; lam = np.array(band) * 1e3
    cen = [np.mean(des.trace(l, c, npup=1500, Csys=Cs)[0]) for l in band]
    rows = []
    for i, l in enumerate(band):
        U, V = des.trace(l, c, npup=2500, Csys=Cs)
        rms = (U - U.mean()).std() * 1e3
        i0, i1 = max(0, i - 1), min(len(band) - 1, i + 1)
        plate = abs((lam[i1] - lam[i0]) / (cen[i1] - cen[i0]))
        rp = d["rB"] / np.cos(beta(d, l) - np.radians(d["beta_ref"]))
        diff = 1.028 * l * (rp / d["EPD"])
        dl = np.hypot(np.hypot(SLIT, 2.355 * rms), np.hypot(PIX, diff)) / 1e3 * plate
        rows.append(dict(lam=int(l * 1e3), plate=float(plate), dl=float(dl), line=float(V.max() - V.min())))
    return rows


def fit_chirp(d, des, Cs):
    """Fit c1,c2,c3 on the FIXED plane (no refocus)."""
    sc = np.array([0.1, 0.01, 1e-4])
    def obj(z):
        out = []
        for l in d["band"]:
            U, _ = des.trace(l, list(z * sc), npup=300, Csys=Cs)
            if len(U) < 50: return 1e6
            out.append((U - U.mean()).std() * 1e3)
        out = np.array(out); return out.mean() + 0.6 * out.max()
    best = min((minimize(obj, s, method="Nelder-Mead",
                         options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
                for s in ([0, 0, 0], [1.8, -1.8, 0], [-1.8, -0.5, 0])), key=lambda r: r.fun)
    return list(best.x * sc)


def fit_c1_only(d, des, Cs):
    from scipy.optimize import minimize_scalar
    def obj(c1):
        out = np.array([(lambda U: (U - U.mean()).std() * 1e3)(des.trace(l, [c1, 0, 0], npup=300, Csys=Cs)[0])
                        for l in d["band"]])
        return out.mean() + 0.6 * out.max()
    return float(minimize_scalar(obj, bounds=(-1.5, 1.5), method="bounded").x)


def idealised_per_lambda(d, des, Cs):
    """Best chirp for EACH wavelength separately (the bound the LPF 'corrected' number assumes)."""
    sc = np.array([0.1, 0.01, 1e-4]); out = []
    for l in d["band"]:
        f = lambda z: (lambda U: (U - U.mean()).std())(des.trace(l, list(z * sc), npup=400, Csys=Cs)[0])
        z = min((minimize(f, s0, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-7, maxiter=1500))
                 for s0 in ([0, 0, 0], [2.5, -1.7, 0], [-2.5, -1.7, 0])), key=lambda r: r.fun).x
        dd = dict(d); dd["band"] = [l]
        out.append(trace_eval_single(d, des, list(z * sc), Cs, l))
    return float(np.mean(out))


def trace_eval_single(d, des, c, Cs, l):
    U, _ = des.trace(l, c, npup=2500, Csys=Cs); rms = (U - U.mean()).std() * 1e3
    b = beta(d, l); r = d["rB"] / np.cos(b - np.radians(d["beta_ref"]))
    plate = (1e6 / d["N"]) * np.cos(b) * np.cos(b - np.radians(d["beta_ref"])) / r   # Shimadzu eq (24)
    return np.hypot(np.hypot(SLIT, 2.355 * rms), np.hypot(PIX, 1.028 * l * r / d["EPD"])) / 1e3 * plate


def run_raytrace(d):
    des = R.Design(d["label"], d["N"], d["alpha"], d["rA"], d["rB"], d["Ry"], d["rho"], d["EPD"], d["band"])
    Cs = plane(d)
    cl = trace_eval(d, des, [0, 0, 0], Cs)
    c1 = fit_c1_only(d, des, Cs)
    one = trace_eval(d, des, [c1, 0, 0], Cs)
    c = fit_chirp(d, des, Cs)
    full = trace_eval(d, des, c, Cs)
    # sign-convention check at mid-band: single-wavelength c1 vs the LPF's defocus-null value
    from scipy.optimize import minimize_scalar
    lm = d["band"][len(d["band"]) // 2]
    c1_mid = minimize_scalar(lambda x: (lambda U: (U - U.mean()).std())(des.trace(lm, [x, 0, 0], npup=400, Csys=Cs)[0]),
                             bounds=(-2, 2), method="bounded").x
    a = np.radians(d["alpha"]); b = beta(d, lm); r = d["rB"] / np.cos(b - np.radians(d["beta_ref"]))
    F20 = np.cos(a)**2/d["rA"] - np.cos(a)/d["Ry"] + np.cos(b)**2/r - np.cos(b)/d["Ry"]
    mean = lambda rows, k: float(np.mean([x[k] for x in rows]))
    return dict(plate=mean(cl, "plate"), astig=mean(cl, "line"), classical=mean(cl, "dl"),
                idealised=idealised_per_lambda(d, des, Cs), one_chirp=mean(one, "dl"), c1=c1,
                full_chirp=mean(full, "dl"), full_min=float(min(x["dl"] for x in full)), cfit=c,
                sign_check=(float(c1_mid), float(-F20 / (lm * 1e-3)), int(lm * 1e3)))


def pct(a, b):
    return abs(a - b) / abs(b) * 100 if b else float("nan")


METRICS = [("plate", "plate scale", "nm/mm"), ("astig", "astig line", "mm"),
           ("classical", "classical res", "nm"), ("idealised", "idealised bound", "nm"),
           ("one_chirp", "one real chirp", "nm")]
NOTES = dict(plate="grating equation", astig="F02 geometry (near-0 for toroid)",
             classical="constant grooves, same fixed plane",
             idealised="defocus nulled per wavelength (not reachable)",
             one_chirp="identical single c1 in both models")


def main():
    rows = []
    for d in DESIGNS:
        rt = run_raytrace(d); lp = run_lpf(d); lp["one_chirp"] = lpf_one_chirp(d, rt["c1"])
        rows.append((d, lp, rt))
    print("=" * 100)
    print("  MODEL CROSS-VALIDATION  |  LPF (analytic) vs ray trace (geometric)  |  identical fixed detector plane")
    print("=" * 100)
    for d, lp, rt in rows:
        print(f"\n{d['label']}   band {int(d['band'][0]*1e3)}-{int(d['band'][-1]*1e3)} nm, f/{d['rA']/d['EPD']:.1f}, "
              f"R={d['Ry']:.2f} mm, plane {d['rB']:.2f} mm @ {d['beta_ref']:+.2f} deg")
        print(f"  {'metric':<15}{'LPF':>11}{'ray trace':>12}{'delta':>10}{'delta %':>10}   note")
        for key, name, unit in METRICS:
            lv, rv = lp[key], rt[key]
            print(f"  {name:<15}{lv:>11.3f}{rv:>12.3f}{rv-lv:>10.3f}{pct(rv,lv):>9.1f}%   {NOTES[key]} ({unit})")
        sc = rt["sign_check"]
        print(f"  one chirp: c1 = {rt['c1']:+.4f} l/mm^2 (both models).  Ray trace with its full c1..c3 fit "
              f"[{rt['cfit'][0]:+.4f}, {rt['cfit'][1]:+.5f}, {rt['cfit'][2]:+.2e}]: mean {rt['full_chirp']:.2f} nm, "
              f"{rt['full_min']:.2f} nm at the defocus null")
        print(f"  sign check @ {sc[2]} nm: trace c1 {sc[0]:+.4f} vs -LPF defocus-null c1 {sc[1]:+.4f}")
    make_figure(rows)
    json.dump([{"design": d["label"], "lpf": lp, "raytrace": rt} for d, lp, rt in rows],
              open("model_validation.json", "w"), indent=1)
    print("\nsaved jasper_model_validation.png and model_validation.json")


def make_figure(rows):
    C_LPF, C_RT = "#1aa5cf", "#e8912a"
    labels = [d["label"].replace("Simulated VLS - ", "Sim ").replace(" (real)", "") for d, _, _ in rows]
    fig, axes = plt.subplots(1, 5, figsize=(19, 4.6))
    fig.suptitle("JASPER model cross-validation — analytic LPF vs geometric ray trace, identical fixed detector plane",
                 fontsize=12.5, fontweight="bold", y=1.03)
    x = np.arange(len(rows)); w = 0.38
    for ax, (key, name, unit) in zip(axes, METRICS):
        lpv = [r[1][key] for r in rows]; rtv = [r[2][key] for r in rows]
        lab_l, lab_r = "LPF (analytic)", "ray trace"
        ax.bar(x - w / 2, lpv, w, label=lab_l, color=C_LPF, edgecolor="#0c5468", linewidth=0.6)
        ax.bar(x + w / 2, rtv, w, label=lab_r, color=C_RT, edgecolor="#8a4d10", linewidth=0.6)
        for xi, (lv, rv) in enumerate(zip(lpv, rtv)):
            ax.text(xi - w / 2, lv, f"{lv:.2f}", ha="center", va="bottom", fontsize=7.5)
            ax.text(xi + w / 2, rv, f"{rv:.2f}", ha="center", va="bottom", fontsize=7.5)
        ax.set_title(f"{name} ({unit})" + ("\nmean over band" if key in ("classical", "idealised", "one_chirp") else "\n"),
                     fontsize=10.5, fontweight="bold")
        ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8, rotation=12, ha="right")
        ax.grid(axis="y", alpha=0.25, linewidth=0.6); ax.set_axisbelow(True)
        top = max(max(lpv), max(rtv)); ax.set_ylim(0, top * 1.45 if top > 0 else 1)
    axes[1].legend(fontsize=8, loc="upper left", frameon=False)
    for k, txt in ((3, "bound: defocus nulled at\nEVERY λ — not reachable\nby any single groove chirp"),
                   (4, "physically real: ONE chirp,\nidentical in both models —\nnulls defocus at ~2 λ only")):
        axes[k].text(0.5, 0.97, txt, transform=axes[k].transAxes, ha="center", va="top", fontsize=7.5,
                     style="italic", color="#8a4d10", bbox=dict(boxstyle="round,pad=0.3", fc="#fff5e9", ec="#e8912a", lw=0.8))
    axes[2].text(0.36, 0.97, "same plane, no free knobs:\ngap = terms the 3rd-order LPF drops\n(higher order; toroid cross-terms)",
                 transform=axes[2].transAxes, ha="center", va="top", fontsize=7.5, style="italic", color="#0c5468",
                 bbox=dict(boxstyle="round,pad=0.3", fc="#e9f7fc", ec="#1aa5cf", lw=0.8))
    fig.tight_layout()
    fig.savefig("jasper_model_validation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
