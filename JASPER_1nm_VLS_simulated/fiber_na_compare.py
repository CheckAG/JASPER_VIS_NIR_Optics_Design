#!/usr/bin/env python3
"""Does a low-NA fiber help? Thorlabs M96L01 (105 um core, NA 0.10) vs Edmund #23-692 (200 um, NA 0.22)
feeding the P0400-01 through a 25 um slit. Traced through the same grating model as relay_design.py.
Relays are idealised here (uniform NA-matched image) so the fibers are compared on equal terms; the real
compact relay reaches ~85 % of its ideal. Fiber power ~ core area x NA^2 (extended diffuse source, e.g. a
reflectance probe on a sample)."""
import numpy as np, relay_design as RD, validate_p0400 as V
rng = np.random.default_rng(5); N = 400_000; SW = 0.025
def cone(NA, n):
    k = NA * np.sqrt(rng.uniform(0, 1, n)); ph = rng.uniform(0, 2 * np.pi, n)
    return np.stack([k * np.cos(ph), k * np.sin(ph), np.sqrt(1 - k * k)], 1)
def disc(R, n):
    r = R * np.sqrt(rng.uniform(0, 1, n)); th = rng.uniform(0, 2 * np.pi, n)
    return np.stack([r * np.cos(th), r * np.sin(th), np.zeros(n)], 1)
def case(core_um, NAf, m, stop):
    """m = relay magnification at the slit (1 = fiber butted to the slit)."""
    Q = disc(core_um * 1e-3 * m / 2, N); D = cone(NAf / m, N)
    keep = np.abs(Q[:, 1]) <= SW / 2; Q, D = Q[keep], D[keep]; slit_frac = keep.mean()
    S, d = RD.to_grating(Q, D); rows = []
    for l in V.BAND:
        U, Vv, ok = RD.grating_trace(S, d, l, stop)
        Vc = np.median(Vv[ok]); on = ok & (np.abs(Vv - Vc) <= 0.1)
        plate, r = RD.plate_scale(l); feff = 10.0 if stop else 1 / (2 * NAf / m)
        fw = np.hypot(RD.fwhm((U[on] - np.median(U[on])) * 1e3), 1.028 * l * feff) / 1e3 * plate
        rows.append((slit_frac * on.sum() / len(U), fw))
    light = np.mean([x[0] for x in rows]) * (core_um / 200) ** 2 * (NAf / 0.22) ** 2
    return light, np.mean([x[1] for x in rows]), max(x[1] for x in rows), min(x[1] for x in rows)
cases = [("#23-692 200 um/0.22 + relay x4.29 (current, ideal relay)", 200, 0.22, 0.22 / 0.0513, False),
         ("M96L01 105 um/0.10 + relay x2 -> f/10",                    105, 0.10, 2.0,  False),
         ("M96L01 105 um/0.10 butted to slit, f/10 grating stop",     105, 0.10, 1.0,  True),
         ("M96L01 105 um/0.10 butted to slit, no stop (f/5)",         105, 0.10, 1.0,  False)]
res = [(n, *case(c, na, m, s)) for n, c, na, m, s in cases]
ref = res[0][1]
print(f"{'configuration':58s} {'light':>7s} {'FWHM mean':>10s} {'range':>13s}")
for n, L, mean, worst, best in res:
    print(f"{n:58s} {L/ref:6.2f}x {mean:8.2f} nm  {best:.2f}-{worst:.2f} nm")
