#!/usr/bin/env python3
"""Standalone sequential ray trace of ../Broadband.OPT — a debugging tool, not a design tool.

Answers one question BeamFour's Layout cannot: *which surface* is killing a ray. When a
wavelength silently vanishes from the layout, run this and it names the blocking surface.

Conventions reproduced from the BeamFour 2.08 sources (B4SOURCES_208):
  - surface frame: pitch P rotates about +Y, so the local axis is (sinP, 0, cosP)
  - the Index column is the medium BEFORE the surface; the medium after surface i is the
    Index of surface i+1
  - grating (RT13.iGrating): in local coords sp = ray.x, and sp_out = sp_in + m*lambda*g,
    with lambda in mm and g in lines/mm. For a +Z ray on a mirror at pitch alpha this gives
    sp_in = -sin(alpha), hence sin(beta) = m*lambda/d - sin(alpha).

Usage:  python3 trace_check.py [n_rays_per_wave]
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE   # release: everything in one directory
OPT = os.path.join(ROOT, "Broadband.OPT")
MED = os.path.join(ROOT, "Broadband.MED")
WAVES = [400, 500, 600, 700, 800, 900, 1000]


def load_med():
    rows = open(MED).read().splitlines()
    hdr = [c.strip() for c in rows[1].split(":")][1:-2]
    lams = [round(float(h) * 1e6) for h in hdr if h]
    idx = {}
    for ln in rows[3:]:
        if not ln.strip():
            continue
        c = [x.strip() for x in ln.split(":")]
        idx[c[0]] = dict(zip(lams, [float(v) for v in c[1:1 + len(lams)]]))
    return idx


def load_opt(path=OPT):
    rows = open(path).read().splitlines()
    out = []
    for ln in rows[3:]:
        if not ln.strip():
            continue
        c = ln.split(":")
        g = lambda i: c[i].strip()
        out.append(dict(mat=g(0), typ=g(1), X=float(g(2)), Z=float(g(3)),
                        pitch=float(g(4)), curv=float(g(5)),
                        gx=float(g(6)) if g(6) else 0.0,
                        order=float(g(7)) if g(7) else 0.0,
                        diam=float(g(8)) if g(8) else 0.0,
                        dx=float(g(9)) if g(9) else 0.0,
                        dy=float(g(10)) if g(10) else 0.0,
                        note=c[12].strip() if len(c) > 12 else ""))
    return out


def to_local(p, d, s):
    ca, sa = math.cos(math.radians(s["pitch"])), math.sin(math.radians(s["pitch"]))
    ex, ez = p[0] - s["X"], p[2] - s["Z"]
    # local x-hat = (cosP,0,-sinP), y-hat = (0,1,0), z-hat = (sinP,0,cosP)
    lp = (ex * ca - ez * sa, p[1], ex * sa + ez * ca)
    ld = (d[0] * ca - d[2] * sa, d[1], d[0] * sa + d[2] * ca)
    return lp, ld


def to_global(lp, ld, s):
    ca, sa = math.cos(math.radians(s["pitch"])), math.sin(math.radians(s["pitch"]))
    x = lp[0] * ca + lp[2] * sa + s["X"]
    z = -lp[0] * sa + lp[2] * ca + s["Z"]
    dx = ld[0] * ca + ld[2] * sa
    dz = -ld[0] * sa + ld[2] * ca
    return (x, lp[1], z), (dx, ld[1], dz)


def intersect(lp, ld, curv, cyl=False):
    """Ray/surface intersection in local coords; returns (point, normal) or None.

    cyl=True treats the surface as a CYLINDER with its axis along local y (curvature
    only in the x-z / dispersion plane) -- a plano-cylindrical field flattener that
    corrects tangential field curvature without touching the sagittal (height) focus.
    """
    if abs(curv) < 1e-12:
        if abs(ld[2]) < 1e-12:
            return None
        t = -lp[2] / ld[2]
        if t < 0:
            return None
        nrm = (0.0, 0.0, 1.0 if ld[2] < 0 else -1.0)
        return tuple(lp[i] + t * ld[i] for i in range(3)), nrm
    R = 1.0 / curv
    if cyl:
        # cylinder axis along local y: solve in x-z only, normal has no y component
        a = ld[0] * ld[0] + ld[2] * ld[2]
        if a < 1e-15:
            return None
        ex, ez = lp[0], lp[2] - R
        b = 2.0 * (ex * ld[0] + ez * ld[2])
        c = ex * ex + ez * ez - R * R
        disc = b * b - 4 * a * c
        if disc < 0:
            return None
        sq = math.sqrt(disc)
        for t in sorted([(-b - sq) / (2 * a), (-b + sq) / (2 * a)], key=abs):
            if t > 1e-9:
                q = tuple(lp[i] + t * ld[i] for i in range(3))
                n = (q[0] / R, 0.0, (q[2] - R) / R)
                nl = math.hypot(n[0], n[2])
                n = (n[0] / nl, 0.0, n[2] / nl)
                if sum(ld[k] * n[k] for k in range(3)) > 0:
                    n = tuple(-v for v in n)
                return q, n
        return None
    e = (lp[0], lp[1], lp[2] - R)
    ed = sum(e[i] * ld[i] for i in range(3))
    disc = ed * ed - (sum(v * v for v in e) - R * R)
    if disc < 0:
        return None
    sq = math.sqrt(disc)
    ts = sorted([-ed - sq, -ed + sq], key=abs)
    for t in ts:
        if t > 1e-9:
            q = tuple(lp[i] + t * ld[i] for i in range(3))
            n = ((q[0]) / R, (q[1]) / R, (q[2] - R) / R)
            nl = math.sqrt(sum(v * v for v in n))
            n = tuple(v / nl for v in n)
            # orient the normal against the incoming ray so cos(incidence) > 0
            if sum(ld[k] * n[k] for k in range(3)) > 0:
                n = tuple(-v for v in n)
            return q, n
    return None


def blocked(s, q):
    r = math.hypot(q[0], q[1])
    if s["dx"] or s["dy"]:                       # rectangular (slit, detector)
        return abs(q[0]) > s["dx"] / 2 or abs(q[1]) > s["dy"] / 2
    if s["diam"]:
        return r > s["diam"] / 2
    return False


def trace(p, d, lam_nm, surfs, med, start=0, want_dir=False):
    """Returns (blocking_surface_index, reason), or (None, local_point) at the detector.

    want_dir=True returns (None, (local_point, local_dir)) so callers can refocus
    analytically -- propagating to a shifted detector plane without re-tracing.
    """
    lam_mm = lam_nm / 1e6
    for i in range(start, len(surfs)):
        s = surfs[i]
        lp, ld = to_local(p, d, s)
        hit = intersect(lp, ld, s["curv"], cyl=("CYL" in s.get("note", "")))
        if hit is None:
            return i, "no intersection"
        q, n = hit
        if blocked(s, q):
            return i, f"aperture (r={math.hypot(q[0], q[1]):.2f})"
        if s["typ"] == "CCD":
            return None, ((q, ld) if want_dir else q)
        if s["typ"] == "iris":
            p, d = to_global(q, ld, s)
            continue
        if s["typ"] == "mirror":
            if s["gx"]:                          # grating, BeamFour iGrating convention
                sp = ld[0] + s["order"] * lam_mm * s["gx"]
                sq_ = ld[1]
                r2 = 1.0 - sp * sp - sq_ * sq_
                if r2 < 0:
                    return i, "RRORD (order beyond 90 deg)"
                ld = (sp, sq_, -math.sqrt(r2))
            else:
                dn = sum(ld[k] * n[k] for k in range(3))
                ld = tuple(ld[k] - 2 * dn * n[k] for k in range(3))
            p, d = to_global(q, ld, s)
            continue
        # lens surface: refract, n_before = this Index, n_after = next surface's Index
        nb = med[s["mat"]][lam_nm] if s["mat"] else 1.0
        nxt = surfs[i + 1] if i + 1 < len(surfs) else None
        na = med[nxt["mat"]][lam_nm] if (nxt and nxt["mat"]) else 1.0
        mu = nb / na
        ci = -sum(ld[k] * n[k] for k in range(3))
        k2 = 1 - mu * mu * (1 - ci * ci)
        if k2 < 0:
            return i, "TIR"
        ld = tuple(mu * ld[k] + (mu * ci - math.sqrt(k2)) * n[k] for k in range(3))
        p, d = to_global(q, ld, s)
    return None, None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    med, surfs = load_med(), load_opt()
    _slit = [i for i, s in enumerate(surfs) if s["note"].startswith("entrance slit")]
    islit = _slit[0] if _slit else -1          # norelay: fiber IS the slit at Z=0, no slit surface
    zslit = surfs[islit]["Z"] if _slit else 0.0

    print(f"{len(surfs)} surfaces; launching from the slit (Z={zslit}) into the accepted cone\n")
    for i, s in enumerate(surfs):
        print(f"  [{i:2}] {s['typ']:6} {s['note'][:62]}")
    print()

    # launch a fan in the dispersion plane, spanning slightly more than the stop accepts
    print(f"{'nm':>5} {'through':>8}  first blocking surface")
    for w in WAVES:
        ok, fails = 0, {}
        for j in range(n):
            u = -0.16 + 0.32 * j / (n - 1)
            d = (u, 0.0, math.sqrt(1 - u * u))
            bad, info = trace((0.0, 0.0, zslit), d, w, surfs, med, start=islit + 1)
            if bad is None:
                ok += 1
            else:
                fails.setdefault(bad, [0, info])
                fails[bad][0] += 1
        if fails:
            worst = max(fails.items(), key=lambda kv: kv[1][0])
            note = surfs[worst[0]]["note"][:44]
            print(f"{w:>5} {ok:>4}/{n}   [{worst[0]}] {note}  x{worst[1][0]} ({worst[1][1]})")
        else:
            print(f"{w:>5} {ok:>4}/{n}   -")


if __name__ == "__main__":
    main()
