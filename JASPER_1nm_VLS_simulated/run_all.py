#!/usr/bin/env python3
"""Run every script in the package in dependency order (see RUN_ORDER.md).

    python run_all.py              # everything (~10-20 min depending on the machine)
    python run_all.py --from 12    # resume from step 12
    python run_all.py --only 12 17 # just these steps

Outputs (PNG / JSON) are written to this folder; the shipped copies live in figures/.
Stops at the first script that fails. Works the same on Windows, macOS and Linux.
"""
import subprocess, sys, time

STEPS = [
    # Release 1 -- simulated custom VLS grating
    (1,  "vls_raytrace.py",           "ray-tracer self-test (console)"),
    (2,  "lpf_resolution.py",         "analytic resolution table -> lpf_all_designs.json"),
    (3,  "lpf_spots.py",              "analytic idealised-bound spots"),
    (4,  "lpf_compare_fig.py",        "what each correction buys (needs 2)"),
    (5,  "raytrace_two_panels.py",    "f/5.3 spherical vs f/8 toroid spots"),
    (6,  "raytrace_three_panels.py",  "f/5.3 -> f/8 -> f/10 + VLS spots (VLS fit)"),
    (7,  "raytrace_layout_light.py",  "layout, light"),
    (8,  "raytrace_layout_scale.py",  "layout, to scale"),
    (9,  "raytrace_layout.py",        "layout, dark"),
    # Release 2 material -- realizability and model cross-check
    (10, "recording_optimize_sim.py", "real 2-point holographic recording -> recording_optimize_sim.json"),
    (11, "validate_models.py",        "analytic vs ray trace on one detector plane"),
    # Release 3 -- Shimadzu P0400-01 practical build
    (12, "validate_p0400.py",         "datasheet reproduction + solved C1-C3 -> validate_p0400.json"),
    (13, "p0400_verify.py",           "pole-vs-edge and orientation tests"),
    (14, "p0400_recording.py",        "lambda_0 back-solve: holographic recordings for the chirp"),
    (15, "p0400_entrance_options.py", "fiber / slit / relay options (needs 12)"),
    (16, "relay_optics.py",           "glass + catalog-lens check (console)"),
    (17, "relay_design.py",           "end-to-end relay design (needs 12) -> relay_design_compact.json"),
    (18, "fiber_na_compare.py",       "0.10 NA fiber comparison (needs 12)"),
    (19, "combined_layout.py",        "full layout, light (needs 12, 17)"),
    (20, "combined_layout_dark.py",   "full layout, dark (needs 12, 17)"),
    # Cross-cutting
    (21, "resolution_metrics.py",     "formula vs measured FWHM / EE80, all designs (needs 10, 12)"),
]

def main(argv):
    first, only = 1, None
    if "--from" in argv: first = int(argv[argv.index("--from") + 1])
    if "--only" in argv: only = {int(x) for x in argv[argv.index("--only") + 1:] if x.isdigit()}
    t_all = time.time()
    for n, script, what in STEPS:
        if n < first or (only and n not in only): continue
        print(f"[{n:2d}] {script:28s} {what}", flush=True); t = time.time()
        r = subprocess.run([sys.executable, script])
        if r.returncode != 0:
            print(f"\nStep {n} ({script}) failed (exit {r.returncode}). Fix it, then: python run_all.py --from {n}")
            sys.exit(r.returncode)
        print(f"     done in {time.time() - t:.0f} s\n", flush=True)
    print(f"All done in {(time.time() - t_all) / 60:.1f} min.")

if __name__ == "__main__":
    main(sys.argv[1:])
