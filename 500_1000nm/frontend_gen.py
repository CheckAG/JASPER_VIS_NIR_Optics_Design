#!/usr/bin/env python3
"""Generate Broadband.OPT — the 500-1000 nm JASPER spectrograph, ready to trace in BeamFour.

Run with no arguments to reproduce the released design exactly:

    python frontend_gen.py            ->  Broadband.OPT

The file it writes is loaded straight into BeamFour (File -> Open). The cylindrical field
flattener is already written in BeamFour's convention (curvature in the Cx column, Curv = 0),
so there is no separate conversion step -- one script, one .OPT file.

Optical layout (global X-Z plane; X = across-dispersion, Z = downstream):

    fiber face = 30 um slit (Z = 0, no relay)
      -> f = 78 mm collimator     (AC254-045-B doublet scaled x1.7333: N-SF6HT flint / N-LAK22 crown)
      -> O16 aperture stop        (f/4.9)
      -> 600 l/mm reflective grating, alpha = 35.64 deg (45 deg deviation mount)
      -> f = 78 mm camera         (same doublet, flipped)
      -> plano-cylindrical N-BK7 field flattener (cyl f = -70 mm, power in the dispersion plane only)
      -> TCD1304 detector, 29.1 mm, at 131 mm along the diffracted arm, tilted -6 deg

Both powered lenses are ONE prescription: the Thorlabs AC254-045-B doublet scaled to f = 78
(curvatures / 1.7333, thicknesses x 1.7333). Scaling preserves the glasses and the shape
factor, so the achromatic correction carries over exactly. The collimator keeps the flint
outer toward the slit; the camera is the same doublet flipped (crown outer to the beam). The
cylinder corrects the camera's dispersion-plane field curvature -- the resolution limiter --
without magnifying the slit-height axis.

Delivered resolution, single exposure: ~1.9 nm mean / 2.9 nm worst across 500-1000 nm, on the
Toshiba TCD1304 (BeamFour-traced; trace_check.py reproduces BeamFour to < 0.1 um per wavelength).

Everything is optional-argument / environment driven for experiments, but the DEFAULTS ARE the
released design:

    python frontend_gen.py [alpha] [stop|-] [det_tilt] [suffix] [flattener_f|-] [det_arm] [relay]

Env overrides (each defaults to the released value): TLO/THI band edges, FCAM camera focal
length, CYLF cylinder focal length, CYLGAP cylinder-to-detector gap, CENTER_NM field decenter.
Passing "relay" restores the AC254-045-B fiber relay in place of the butted fiber.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE   # release: everything in one directory

# ---- CLI parameters ---------------------------------------------------------
# RELEASE DEFAULTS reproduce the validated 500-1000 nm / 1.9 nm design with no arguments:
#   fiber butted to the slit (no relay), f=78 collimator+camera, O16 stop, 600 l/mm grating,
#   cylindrical field flattener (f=-70), detector at 131 mm along the arm, tilt -6 deg.
# Pass "relay" to restore the AC254-045-B fiber relay for experiments.
NORELAY = "relay" not in sys.argv
sys.argv = [a for a in sys.argv if a not in ("norelay", "relay")]
ALPHA = float(sys.argv[1]) if len(sys.argv) > 1 else 35.64
STOP = None if (len(sys.argv) < 3 or sys.argv[2] in ("-", "none", "0")) else float(sys.argv[2])
DET_TILT = float(sys.argv[3]) if len(sys.argv) > 3 else -6.0
SUFFIX = sys.argv[4] if len(sys.argv) > 4 else ""
# 5th arg: field-flattener focal length in mm (negative). "-" or absent = no flattener.
FLAT_F = None if (len(sys.argv) < 6 or sys.argv[5] in ("-", "none", "0")) else float(sys.argv[5])
# 6th arg: detector distance along the arm from the grating. Each camera configuration has its
# own best focus, so a flattener variant MUST re-specify this or it will be badly defocused.
DET_ARM_CLI = float(sys.argv[6]) if len(sys.argv) > 6 else None
OPTFILE = os.path.join(ROOT, f"Broadband{SUFFIX}.OPT")

if len(sys.argv) < 3:
    STOP = 16.0

# ---- design constants -------------------------------------------------------
LPMM = 600
# The camera and detector are centred on the ray that bisects the SPECTRUM, not on the
# middle wavelength. Because beta(lambda) is nonlinear, 700 nm sits 0.2 mm off the midpoint
# of the 400..1000 landing points. Centring on beta_c = (beta400 + beta1000)/2 equalises the
# detector margins instead of leaving 1.35 / 1.76 mm.
LAM_LO = float(os.environ.get("TLO", "500"))    # band (env-tunable for grating-rotated sub-bands)
LAM_HI = float(os.environ.get("THI", "1000"))
SCALE = float(os.environ.get("FCAM", "78")) / 45.0   # camera focal length (env-tunable for the flattener design)
# The two arms need different apertures. The collimator only ever sees the O16 beam plus the
# 1 mm slit height (~17 mm), so O25.4 is ample. Only the CAMERA needs O40: it sits in the
# dispersed fan, where the footprint is 35.6 mm.
DIAM_COLL = 25.4                    # collimator (house convention for a O25.4 barrel is 22.86 clear)
DIAM_CAM = 38.0                     # O40 barrel, 38 mm clear
DIAM_GRT = 50.0                     # 50x50 mm grating; footprint is 24.6 mm
DET_DX, DET_DY = 29.100, 0.200      # TCD1304 3648 x 8 um
SLIT_DX, SLIT_DY = 0.030, 1.000     # 30 um x 1 mm entrance slit

D_SLIT_COLL = 40.6037 * SCALE       # slit -> collimator s1 (scaled AC254-045-B BFL)
D_COLL_GRT = 70.000                 # collimator s3 -> grating
D_STOP = 20.000                     # collimator s3 -> aperture stop
# Grating -> camera s1. Squeezed from both sides:
#   too near  -> the O38 lens body cuts into the incoming O16 collimated beam
#   too far   -> the dispersed fan overruns the 38 mm clear aperture
# At the O16 stop, 44.5 mm gives +4.6 mm and +2.4 mm against the two walls. Note the stop
# had to come down from O20 to O16: at O20 the fan footprint is 40.4 mm and overruns the
# lens, which is what clipped 700-1000 nm out of the first trace.
D_GRT_CAM = 44.500
# Detector distance measured along the arm FROM THE GRATING (not from the last lens), because
# that is the focus parameter. 131.0 mm with tilt -6.0 deg is the released optimum for the
# f=78 camera + cylindrical flattener: ~1.9 nm mean / 2.9 nm worst, 28.37 mm spectrum extent,
# no vignetting. A different camera or flattener has its own best focus, so any variant must
# re-specify this (6th CLI arg) or it will be defocused.
DET_ARM = 131.000                   # released optimum for the f=78 cylinder design
N_BK7_700 = 1.51306                 # flattener glass, index at band centre
FLAT_GAP = 8.333                    # flattener s1 -> detector, mm (optimised)
FLAT_T = 2.500                      # flattener centre thickness, mm
if DET_ARM_CLI:
    DET_ARM = DET_ARM_CLI

T_FLINT = 1.8 * SCALE               # 3.120 at f=78
T_CROWN = 6.0 * SCALE               # 10.400 at f=78

C1 = 0.007859 / SCALE               # weak flint outer   (R  127.2 -> 197.9)
C2 = 0.039920 / SCALE               # cemented interface (R   25.1 ->  39.0)
C3 = 0.034037 / SCALE               # strong crown outer (R  -29.4 -> -45.7)

FRONT = [
    ("", "lens", 0.0, 40.6037, 0.0, 0.007859, "", "", 22.860, "", "", "", "relay L1 AC254-045-B s1 (N-SF6HT flint outer)"),
    ("N-SF6HT", "lens", 0.0, 42.4037, 0.0, 0.039920, "", "", 22.860, "", "", "", "relay L1 s2 (N-SF6HT -> N-LAK22)"),
    ("N-LAK22", "lens", 0.0, 48.4037, 0.0, -0.034037, "", "", 22.860, "", "", "", "relay L1 s3 (N-LAK22 -> gap/filter slot)"),
    ("", "lens", 0.0, 58.4037, 0.0, 0.034037, "", "", 22.860, "", "", "", "relay L2 AC254-045-B s1 (N-LAK22 crown outer)"),
    ("N-LAK22", "lens", 0.0, 64.4037, 0.0, -0.039920, "", "", 22.860, "", "", "", "relay L2 s2 (N-LAK22 -> N-SF6HT)"),
    ("N-SF6HT", "lens", 0.0, 66.2037, 0.0, -0.007859, "", "", 22.860, "", "", "", "relay L2 s3 (N-SF6HT -> slit)"),
    ("", "iris", 0.0, 106.8074, 0.0, 0.0, "", "", "", SLIT_DX, SLIT_DY, "S",
     f"entrance slit {SLIT_DX * 1000:.0f}um x {SLIT_DY:.0f}mm"),
]
Z_SLIT = 106.8074

if NORELAY:
    # The fiber face IS the slit, at the origin. The slit is the ray SOURCE, not an obstruction,
    # so it is not emitted as a surface -- a ray starting in the plane of an iris is a
    # degenerate (t=0) intersection. Its 24.7 um image stays in the budget analytically, the
    # same way build_spot_html.py already accounts for it.
    FRONT = []
    Z_SLIT = 0.0


def row(mat, typ, X, Z, pitch, curv, gx, order, diam, dx, dy, form, note):
    """One optics-table line, including BeamFour's Cx column.

    A cylindrical surface (its note contains CYL, and it is the curved side, not the plano
    side) is written the way BeamFour expects a cylinder: its curvature goes in the Cx column
    with Curv = 0, so it is curved in the x-z dispersion plane and flat along the slit height.
    Every other surface leaves Cx blank and keeps its curvature in Curv. Emitting Cx here,
    rather than as a separate post-processing pass, keeps the pipeline at a single .OPT file."""
    def num(v, w, p):
        return (" " * w) if v == "" else f"{v:{w}.{p}f}"
    if "CYL" in note and "plano" not in note:
        cx = f"{curv:10.6f}"        # curvature -> Cx (power in the dispersion plane only)
        curv = 0.0                  # Curv = 0  (flat along the slit height)
    else:
        cx = " " * 10               # blank Cx: ordinary spherical or flat surface
    return (
        f"{mat:>10}:{typ:>6}:{X:9.4f}:{Z:9.4f}:{pitch:9.4f}:{curv:10.6f}:{cx}:"
        f"{str(gx):>5}:{str(order):>6}:{num(diam,7,3)}:{num(dx,7,3)}:{num(dy,7,3)}:"
        f"{form:>4}: {note}"
    )


def build():
    rows = list(FRONT)

    # collimator: scaled AC254-045-B, flint outer facing the slit (relay-L1 orientation)
    zc1 = Z_SLIT + D_SLIT_COLL
    zc2 = zc1 + T_FLINT
    zc3 = zc2 + T_CROWN
    rows += [
        ("", "lens", 0.0, zc1, 0.0, C1, "", "", DIAM_COLL, "", "", "",
         "collimator AC254-045-B x1.7333 s1 (N-SF6HT flint outer, f=78 O25.4)"),
        ("N-SF6HT", "lens", 0.0, zc2, 0.0, C2, "", "", DIAM_COLL, "", "", "",
         "collimator s2 (N-SF6HT -> N-LAK22)"),
        ("N-LAK22", "lens", 0.0, zc3, 0.0, -C3, "", "", DIAM_COLL, "", "", "",
         "collimator s3 (N-LAK22 -> collimated space)"),
    ]

    # aperture stop in the collimated space: sets both f/# and the camera field footprint
    if STOP:
        rows.append(("", "iris", 0.0, zc3 + D_STOP, 0.0, 0.0, "", "", STOP, "", "", "",
                     f"aperture stop O{STOP:.0f} (camera f/{SCALE * 45.0 / STOP:.1f}, collimator f/{SCALE * 45.0 / STOP:.1f})"))

    # grating (reflective) on the +Z axis
    zg = zc3 + D_COLL_GRT
    rows.append(("", "mirror", 0.0, zg, ALPHA, 0.0, LPMM, 1, DIAM_GRT, "", "", "S",
                 f"grating {LPMM} l/mm order1 REFLECTIVE, alpha={ALPHA} (45 deg deviation mount)"))

    # the arm follows the ray that bisects the spectrum (see LAM_LO/LAM_HI above)
    d_nm = 1.0e6 / LPMM
    sa = math.sin(math.radians(ALPHA))

    def bet(lam):
        return math.degrees(math.asin(lam / d_nm - sa))

    _center = os.environ.get("CENTER_NM")   # decenter the field: put this wavelength on-axis (monotonic focus)
    beta = bet(float(_center)) if _center else (bet(LAM_LO) + bet(LAM_HI)) / 2.0
    lam_c = d_nm * (math.sin(math.radians(beta)) + sa)
    theta = -(180.0 - (ALPHA - beta))
    ux, uz = math.sin(math.radians(theta)), math.cos(math.radians(theta))

    # the two walls the camera distance sits between (see D_GRT_CAM above).
    #
    # The footprint must use the DIFFRACTED beam width, not the incident one. A reflective
    # grating is anamorphic: it illuminates a length beam/cos(alpha) and re-emits it as
    # beam*cos(beta)/cos(alpha), so at alpha=35.64 the O16 beam leaves 19-20 mm wide. Sizing
    # the lens on the incident diameter under-reads the footprint by ~4 mm and silently
    # clips the band edges -- which is exactly what killed 700-1000 nm in the first trace.
    r_beam = (STOP / 2.0) if STOP else DIAM_CAM / 2.0
    x_cam = D_GRT_CAM * ux
    clearance = abs(x_cam) - r_beam - DIAM_CAM / 2
    edges = []
    for lam in (LAM_LO, 500.0, 600.0, 700.0, 800.0, 900.0, LAM_HI):
        bl = bet(lam)
        off = D_GRT_CAM * math.tan(math.radians(bl - beta))
        hw = r_beam * math.cos(math.radians(bl)) / math.cos(math.radians(ALPHA))
        edges += [off - hw, off + hw]
    footprint = max(edges) - min(edges)

    print(f"# alpha={ALPHA} beta_c={beta:.3f} ({lam_c:.1f} nm) deviation={ALPHA - beta:.3f} arm={theta:.3f}")
    print(f"# camera s1 at X={x_cam:.1f}mm, {D_GRT_CAM} mm from grating")
    print(f"#   beam clearance  {clearance:+.1f} mm  (lens body vs the O{2 * r_beam:.0f} incoming beam)")
    print(f"#   fan footprint   {footprint:.1f} mm  ({DIAM_CAM - footprint:+.1f} mm of {DIAM_CAM:.0f} clear,"
          f" diffracted width included)")
    if clearance < 0 or footprint > DIAM_CAM:
        print("#   *** camera arm does not fit -- shrink the stop, move D_GRT_CAM, or open the lens ***")

    def along(dist):
        return dist * ux, zg + dist * uz

    # camera: same scaled doublet, flipped (crown outer to the collimated beam)
    specs = [
        (D_GRT_CAM, C3,  "",        "camera AC254-045-B x1.7333 s1 (N-LAK22 crown outer, f=78 O40)"),
        (T_CROWN,   -C2, "N-LAK22", "camera s2 (N-LAK22 -> N-SF6HT)"),
        (T_FLINT,   -C1, "N-SF6HT", "camera s3 (N-SF6HT -> detector)"),
    ]
    dist = 0.0
    for step, curv, mat, note in specs:
        dist += step
        X, Z = along(dist)
        rows.append((mat, "lens", X, Z, theta, curv, "", "", DIAM_CAM, "", "", "", note))

    # optional field flattener: plano-concave, concave toward the camera, plano toward the
    # detector. Lensmaker with a plano back gives c1 = 1/(f*(n-1)); f<0 -> c1<0.
    if FLAT_F:
        c_flat = 1.0 / (FLAT_F * (N_BK7_700 - 1.0))
        X, Z = along(DET_ARM - FLAT_GAP)
        rows.append(("", "lens", X, Z, theta, c_flat, "", "", DIAM_CAM, "", "", "",
                     f"field flattener s1 (PCV f={FLAT_F:.0f}, R={1/c_flat:.1f}, concave->camera)"))
        X, Z = along(DET_ARM - FLAT_GAP + FLAT_T)
        rows.append(("N-BK7", "lens", X, Z, theta, 0.0, "", "", DIAM_CAM, "", "", "",
                     "field flattener s2 (plano -> detector)"))

    # CYLINDRICAL field flattener (env-driven, ON by default): power only in the dispersion (x)
    # plane -- flattens the tangential focal surface, the only one that sets resolution, without
    # magnifying the height. row() writes it in BeamFour's cylinder convention (curvature in the
    # Cx column, Curv = 0); trace_check.py reads that same Cx column. Its curvature is carried in
    # the normal `curv` slot here and moved to Cx at write time. CYLF = cyl focal length mm (neg).
    _cylf = os.environ.get("CYLF", "-70")   # cylinder ON by default
    if _cylf:
        _cylf = float(_cylf)
        _cgap = float(os.environ.get("CYLGAP", "15"))
        c_cyl = 1.0 / (_cylf * (N_BK7_700 - 1.0))
        X, Z = along(DET_ARM - _cgap)
        rows.append(("", "lens", X, Z, theta, c_cyl, "", "", DIAM_CAM, "", "", "",
                     f"field flattener CYL s1 (cyl f={_cylf:.0f}, R={1/c_cyl:.1f}, concave->camera)"))
        X, Z = along(DET_ARM - _cgap + FLAT_T)
        rows.append(("N-BK7", "lens", X, Z, theta, 0.0, "", "", DIAM_CAM, "", "", "",
                     "field flattener CYL s2 (plano -> detector)"))

    X, Z = along(DET_ARM)
    rows.append(("", "CCD", X, Z, theta + DET_TILT, 0.0, "", "", "", DET_DX, DET_DY, "S",
                 f"detector TCD1304 {DET_DX}mm x {DET_DY}mm, {DET_ARM} mm along arm, tilt {DET_TILT}"))
    return rows


def main():
    rows = build()
    n = len(rows)
    tag = f"alpha={ALPHA}" + (f", stop O{STOP:.0f}" if STOP else ", no stop") + f", tilt {DET_TILT}"
    front = "fiber butted to 30um slit at Z=0" if NORELAY else "AC254-045-B relay + 30um slit"
    header = [
        f"{n}  {os.path.basename(OPTFILE)}  500-1000 nm spectrograph ({tag}): {front} "
        f"-> f78 doublet / 600 l/mm REFLECTIVE / f78 doublet / cyl flattener -> TCD1304",
        "     Index:  type:        X:        Z:    Pitch:      Curv:        Cx:   Gx: Order:   Diam:     Dx:     Dy:Form: notes",
        "----------:------:---------:---------:---------:----------:----------:-----:------:-------:-------:-------:----:",
    ]
    with open(OPTFILE, "w") as f:
        f.write("\n".join(header + [row(*r) for r in rows]) + "\n")
    print(f"wrote {OPTFILE}  ({n} surfaces)")


if __name__ == "__main__":
    main()
