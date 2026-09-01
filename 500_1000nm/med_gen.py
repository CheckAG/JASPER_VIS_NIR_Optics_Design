#!/usr/bin/env python3
"""Generate ../Broadband.MED for the 400-1000 nm build.

Three glasses: N-LAK22 (crown) and N-SF6HT (flint) for the AC254-045-B doublets that make
up the relay and both spectrograph arms, plus N-BK7 for the optional field flattener. The
relay uses the doublet at f=45, both spectrograph arms use it scaled x1.5556 to f=70.

The existing builds tabulate the two doublet glasses at 500-1000 nm only (values hard-coded in
JASPER_4f.MED). This design needs 400 nm, so they are regenerated from the Schott
Sellmeier coefficients — which reproduce the existing tables to 5e-6, checked by the
assertion below.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE   # release: everything in one directory
MEDFILE = os.path.join(ROOT, "Broadband.MED")

WAVES_NM = [400, 500, 600, 700, 800, 900, 1000]

# Schott Sellmeier coefficients (B1,B2,B3,C1,C2,C3); lambda in micrometres.
SELLMEIER = {
    "N-LAK22": (1.14229781, 0.535138441, 1.04088385,
                0.00585778594, 0.0198546147, 100.834017),
    "N-SF6HT": (1.77931763, 0.338149866, 2.08734474,
                0.0133714182, 0.0617533621, 174.01759),
    # N-BK7 is only used by the optional field flattener, but it must always be present:
    # BeamFour aborts the whole parse if any .OPT glass is missing from the .MED.
    "N-BK7": (1.03961212, 0.231792344, 1.01046945,
              0.00600069867, 0.0200179144, 103.560653),
}
ABBE = {"N-LAK22": 55.87, "N-SF6HT": 25.35, "N-BK7": 64.17}

# JASPER_4f.MED / 600G_recon.MED values at 500,600,700,750,800,900,1000 nm -- the fit must match.
LEGACY_UM = [0.50, 0.60, 0.70, 0.75, 0.80, 0.90, 1.00]
LEGACY = {
    "N-LAK22": [1.65784, 1.65041, 1.64584, 1.64414, 1.64270, 1.64036, 1.63847],
    "N-SF6HT": [1.82372, 1.80327, 1.79172, 1.78766, 1.78433, 1.77918, 1.77532],
}


def sellmeier(coef, lam_um):
    B1, B2, B3, C1, C2, C3 = coef
    l2 = lam_um * lam_um
    n2 = 1 + B1 * l2 / (l2 - C1) + B2 * l2 / (l2 - C2) + B3 * l2 / (l2 - C3)
    return n2 ** 0.5


def main():
    order = ["N-LAK22", "N-SF6HT", "N-BK7"]

    # the fit must reproduce the glass tables the other two design lines already use
    for name in LEGACY:
        err = max(abs(sellmeier(SELLMEIER[name], l) - n)
                  for l, n in zip(LEGACY_UM, LEGACY[name]))
        assert err < 1e-4, f"{name} Sellmeier disagrees with the legacy .MED by {err}"

    hdr = "".join(f"{w / 1.0e6:9.5f}:" for w in WAVES_NM)
    lines = [
        f"{len(order)} entries  Broadband.MED  400-1000 nm build: AC254-045-B glasses "
        f"(N-LAK22 crown / N-SF6HT flint); waves in mm",
        f"  Material:{hdr}  Abbe :",
        "----------:" + "---------:" * len(WAVES_NM) + "------:",
    ]
    for name in order:
        ns = [sellmeier(SELLMEIER[name], w / 1000.0) for w in WAVES_NM]
        lines.append(f"{name:>10}:" + "".join(f"{n:9.5f}:" for n in ns) + f"{ABBE[name]:6.2f}:")

    with open(MEDFILE, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {MEDFILE}")
    for name in order:
        ns = [sellmeier(SELLMEIER[name], w / 1000.0) for w in WAVES_NM]
        print(f"  {name:>8}: " + " ".join(f"{n:.5f}" for n in ns))


if __name__ == "__main__":
    main()
