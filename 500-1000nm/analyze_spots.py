#!/usr/bin/env python3
"""Parse the BeamFour-traced Broadband.RAY into per-(wavelength, field) spot statistics.

Workflow: in BeamFour, load Broadband_B4.OPT + Broadband.MED + Broadband.RAY, Run -> InOut,
then SAVE the ray table back to Broadband.RAY *in this directory*. This script reads it and
writes spots.json (consumed by build_spot_html.py) next to itself.

It groups the surviving rays by (wavelength, field) and computes the centroid, RMS radius
(x = dispersion axis, y = slit-height axis, r = radial) and the geometric (max) radius, in um.

Column layout of the InOut-saved .RAY (fixed width, colon-delimited):
    X0 Y0 U0 V0 @wave ygoal xfinal yfinal notes
    wave tag at chars 40-51, xfinal at 60-70, yfinal at 71-81, field label is the last token
    of the note. A blank xfinal means the ray was vignetted/lost (skipped).
"""
import collections, math, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
RAYFILE = os.path.join(HERE, "Broadband.RAY")   # release: same directory as this script
OUTFILE = os.path.join(HERE, "spots.json")

NM = {"0.00050": 500, "0.00060": 600, "0.00070": 700, "0.00075": 750,
      "0.00080": 800, "0.00090": 900, "0.00100": 1000}
WAVES = [500, 600, 700, 750, 800, 900, 1000]
FIELDS = ["ctr", "+edge", "-edge"]


def parse():
    print("reading:", RAYFILE)
    rows = open(RAYFILE).read().splitlines()[3:]
    hits = collections.defaultdict(list)      # (nm, field) -> [(x_mm, y_mm)]
    launched = collections.defaultdict(int)
    for ln in rows:
        if len(ln) < 52:
            continue
        w = NM.get(ln[40:51].strip())
        note = ln[71:].strip()
        field = note.split()[-1] if note else "?"
        launched[(w, field)] += 1
        try:
            hits[(w, field)].append((float(ln[60:70]), float(ln[71:81])))
        except ValueError:
            pass                              # blank xfinal = vignetted/lost
    return hits, launched


def stats(pts):
    n = len(pts)
    cx = sum(p[0] for p in pts) / n
    cy = sum(p[1] for p in pts) / n
    rx = math.sqrt(sum((p[0] - cx) ** 2 for p in pts) / n) * 1000
    ry = math.sqrt(sum((p[1] - cy) ** 2 for p in pts) / n) * 1000
    rr = math.hypot(rx, ry)
    geo = max(math.hypot(p[0] - cx, p[1] - cy) for p in pts) * 1000
    return dict(cx=round(cx, 4), cy=round(cy, 4), rmsx=round(rx, 1), rmsy=round(ry, 1),
                rmsr=round(rr, 1), geo=round(geo, 1),
                pts=[[round((p[0] - cx) * 1000, 2), round((p[1] - cy) * 1000, 2)] for p in pts])


def main():
    hits, launched = parse()
    out = {}
    hdr = f"{'nm':>5} {'fld':>6} {'n':>8} {'cx(mm)':>9} {'RMSx(um)':>9} {'GEO(um)':>8}"
    print(hdr)
    for w in WAVES:
        for f in FIELDS:
            L = launched[(w, f)]; pts = hits[(w, f)]; n = len(pts)
            if n == 0:
                out[f"{w}_{f}"] = dict(n=0, L=L)
                if f == "ctr":
                    print(f"{w:>5} {f:>6} {0:>3}/{L:<4}")
                continue
            st = stats(pts); st.update(n=n, L=L)
            out[f"{w}_{f}"] = st
            if f == "ctr":
                print(f"{w:>5} {f:>6} {n:>3}/{L:<4} {st['cx']:>9.3f} {st['rmsx']:>9.1f} {st['geo']:>8.1f}")
    json.dump(out, open(OUTFILE, "w"))
    span = out["500_ctr"]["cx"] - out["1000_ctr"]["cx"]
    tag = "CORRECT 131mm design" if abs(span - 28.37) < 0.15 else "*** span off -- retraced the wrong OPT? ***"
    print(f"\nspan {span:.2f} mm  ->  {tag}")
    print(f"wrote {OUTFILE}")


if __name__ == "__main__":
    main()
