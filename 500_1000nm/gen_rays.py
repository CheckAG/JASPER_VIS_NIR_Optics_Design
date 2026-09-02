#!/usr/bin/env python3
"""Generate Broadband.RAY for the 500-1000 nm build (fiber-launch, colour-tagged).

    python gen_rays.py fan            # 5-ray fan/wavelength (chief + 4 marginals) -> Layout
    python gen_rays.py spots [rings]  # 3 fields x hexapolar pupil x waves -> spot diagrams

Rays launch from the fiber face at the slit plane (Z=0) into the NA 0.22 emission cone and
enter the 500-1000 nm 600 l/mm reflective spectrograph directly (no relay in the released
design). Overwrites Broadband.RAY in this directory.

*** HARD LIMIT: BeamFour reads at most JMAX = 3600 rows. ***

B4constants.java:59 sets `JMAX = 3600` for the ray table, and REJIF sizes wavenames[] from it.
A longer .RAY is NOT rejected -- it is silently TRUNCATED. Because rays are emitted grouped by
wavelength, an overlong file drops whole wavelengths off the end of the list: a 9849-ray file
loaded only 400 nm, 500 nm and part of 600 nm, and 700-1000 nm vanished from the layout with no
error message. That cost a long debugging detour; do not raise MAX_ROWS to "get more rays".

So the ray budget is fixed and the only question is how to spend it. Two levers:

  * pupil rings   -> 1 + 6*(1+2+...+k) points, the thing that actually sets spot fidelity
  * field points  -> must lie INSIDE the detector, see FIELDS_MM below

Rings vs cost, at 7 waves x 3 fields (21 groups):

    rings   pupil pts   total rays   fits in 3600?
      6         127         2667      yes
      7         169         3549      yes  <- the most this layout can afford
      8         217         4557      NO   (silently drops 700-1000 nm)
     12         469         9849      NO

For a denser pupil, cut fields rather than exceeding the limit: 1 field x 7 waves leaves room
for rings=12 (469 pts, 3283 rays). generate() enforces the cap and says what it did.
"""
import sys
import os
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))   # release: one directory
RAYFILE = os.path.join(ROOT, "Broadband.RAY")

MAX_ROWS = 3600     # BeamFour B4constants.JMAX -- files longer than this are silently truncated
HEADER_ROWS = 3     # title + column header + ruler, all counted against JMAX

COLOR = {"0.00050": "b", "0.00060": "c", "0.00070": "g",
         "0.00080": "y", "0.00090": "m", "0.00100": "r"}
WAVES = ["0.00050", "0.00060", "0.00070", "0.00080", "0.00090", "0.00100"]  # 500-1000 nm build
NA = 0.12           # 100 um low-NA fiber emission cone (nearly matches the f/4.9 stop)
RINGS = 7           # the most a 7-wave x 3-field set fits under JMAX; see the docstring
# Measured, not assumed: binary-searched the largest NA at the slit that reaches the detector
# (trace_check), 0.1136 at 400 nm .. 0.1145 at 700 nm -- i.e. exactly the f/4.4 the O16 stop
# sets. Sampling the full 0.22 fiber cone throws two thirds of a fixed 3600-row budget at rays
# the stop kills. 0.125 keeps a visible vignetting margin and spends the rest on real spots.
NA_SPOTS = 0.11     # accepted cone at the f=78 camera through the O16 stop (~f/4.9)
NA_STOP = 0.102

# Field points along the slit, in mm at the fiber. The detector is only 0.200 mm TALL, and the
# slit->detector magnification in y is 1.0, so the old +-0.5 mm ("slit edge") fields landed
# +-0.5 mm off a +-0.1 mm detector and could never register. +-0.08 mm samples the slit height
# that the TCD1304 actually sees.
FIELDS_MM = [(0.0, "ctr"), (0.05, "+edge"), (-0.05, "-edge")]  # 100 um fiber core, within the 0.2 mm detector height

RHEADER = "   X0    :   Y0    :   U0    :   V0    :  @wave    : ygoal :  xfinal  :  yfinal  : notes"
RRULER = "---------:---------:---------:---------:-----------:-------:----------:----------:"


def line(x0, y0, u0, v0, wave, note):
    wf = f"{float(wave):11.5f}" + COLOR[wave]
    return (f"{x0:9.3f}:{y0:9.3f}:{u0:9.4f}:{v0:9.4f}:{wf}"
            f"{0.0:7.3f}:{'':>10}:{'':>10}: {note}")


def hexapolar(na=NA, rings=6):
    pts = [(0.0, 0.0)]
    for k in range(1, rings + 1):
        r = na * k / rings
        for j in range(6 * k):
            a = 2 * np.pi * j / (6 * k)
            pts.append((r * np.cos(a), r * np.sin(a)))
    return pts


def build_fan():
    dirs = [(0.0, 0.0, "chief"), (0.22, 0.0, "+U marg"), (-0.22, 0.0, "-U marg"),
            (0.0, 0.22, "+V marg"), (0.0, -0.22, "-V marg")]
    lines = []
    for w in WAVES:
        nm = int(round(float(w) * 1e6))
        for u0, v0, lbl in dirs:
            lines.append(line(0.0, 0.0, u0, v0, w, f"{nm}nm {lbl}"))
    return lines, "5-ray fan (chief + 4 marginals @ NA0.22) x 7 waves"


def max_rings(nfields, nwaves=len(WAVES)):
    """Largest ring count whose ray total still fits under BeamFour's JMAX."""
    budget = MAX_ROWS - HEADER_ROWS
    k = 1
    while len(hexapolar(rings=k + 1)) * nfields * nwaves <= budget:
        k += 1
    return k


def build_spots(rings=RINGS):
    fields = FIELDS_MM
    pts = hexapolar(na=NA_SPOTS, rings=rings)
    lines = []
    for w in WAVES:
        nm = int(round(float(w) * 1e6))
        for y0, fl in fields:
            for u0, v0 in pts:
                lines.append(line(0.0, y0, u0, v0, w, f"{nm} {fl}"))
    return lines, (f"spot-diagram set: {len(fields)} fields x {len(pts)} pupil pts "
                   f"({rings}-ring hexapolar NA{NA_SPOTS}) x {len(WAVES)} waves")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "spots"
    rings = int(sys.argv[2]) if len(sys.argv) > 2 else RINGS
    cap = max_rings(len(FIELDS_MM))
    if mode != "fan" and rings > cap:
        need = len(hexapolar(rings=rings)) * len(FIELDS_MM) * len(WAVES)
        print(f"  ! rings={rings} needs {need} rays; BeamFour reads only {MAX_ROWS - HEADER_ROWS}"
              f" and SILENTLY TRUNCATES, dropping whole wavelengths.")
        print(f"  ! capped to rings={cap}. For a denser pupil, use fewer fields.")
        rings = cap
    lines, desc = build_fan() if mode == "fan" else build_spots(rings)
    title = f"{len(lines)}  Broadband.RAY  {desc}"
    rows = HEADER_ROWS + len(lines)
    assert rows <= MAX_ROWS, (
        f"{rows} rows exceeds BeamFour's JMAX={MAX_ROWS}; it would truncate silently and "
        f"drop whole wavelengths from the trace")
    with open(RAYFILE, "w") as f:
        f.write(title + "\n" + RHEADER + "\n" + RRULER + "\n" + "\n".join(lines) + "\n")
    print(f"[{mode}] wrote {len(lines)} rays ({rows}/{MAX_ROWS} rows) -> {RAYFILE}")
    if mode != "fan":
        pup = hexapolar(na=NA_SPOTS, rings=rings)
        keep = sum(1 for u, v in pup if (u * u + v * v) ** 0.5 <= NA_STOP)
        print(f"  {len(pup)} pupil pts/field-wave; ~{keep} inside the O16 stop"
              f" ({100 * keep / len(pup):.0f}%) -> ~{keep} pts per spot")
        print(f"  all {len(WAVES)} wavelengths fit: {len(lines) // len(WAVES)} rays each")


if __name__ == "__main__":
    main()
