#!/usr/bin/env python3
"""Convert Broadband.OPT -> Broadband_B4.OPT (BeamFour cylindrical convention).

Our frontend_gen.py writes the field flattener with a "CYL" note and its curvature in the
normal Curv column (which trace_check.py reads directly). BeamFour instead expects a cylinder
as Curv = 0 with the curvature in a **Cx** column (curved in x-z, flat in y). This script moves
the flattener's curvature from Curv into a new Cx column and sets Curv = 0, producing the file
you actually load into BeamFour. Everything stays in this directory.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "Broadband.OPT")
DST = os.path.join(HERE, "Broadband_B4.OPT")


def insert_after_curv(line, field):
    parts = line.split(':')
    return ':'.join(parts[:6] + [field] + parts[6:])   # Curv is column index 5


def main():
    src = open(SRC).read().splitlines()
    out = [src[0],                                   # title
           insert_after_curv(src[1], '        Cx'),  # header
           insert_after_curv(src[2], '----------')]  # ruler
    for ln in src[3:]:
        if not ln.strip():
            continue
        p = ln.split(':')
        note = p[12] if len(p) > 12 else ''
        if 'CYL' in note and 'plano' not in note:   # the curved cylinder surface
            cx = f'{float(p[5]):10.6f}'
            p[5] = f'{0.0:10.6f}'                    # Curv = 0 (flat in y)
        else:
            cx = '          '                        # blank Cx
        out.append(':'.join(p[:6] + [cx] + p[6:]))
    open(DST, "w").write('\n'.join(out) + '\n')
    print(f"wrote {DST}")


if __name__ == "__main__":
    main()
