# JASPER — 500–1000 nm VIS-NIR spectrometer (BeamFour design)

A fully reproducible optical design for a low-cost, right-to-repair fiber-fed grating
spectrometer covering **500–1000 nm**, built from stock 1-inch optics plus one custom
cylindrical corrector, and validated in **BeamFour** (a free geometric ray tracer).

**Delivered resolution (BeamFour-traced, single exposure): ~1.9 nm mean, 2.9 nm worst**
across 500–1000 nm on a Toshiba TCD1304 (29.1 mm) detector. What it would take to push
further, toward **~1 nm**, is described below.

Everything in this folder is self-contained: the Python scripts generate the BeamFour input
files, you trace them in BeamFour, and the scripts turn the trace into a spot-analysis report.
**All files live in this one directory** — scripts read and write here, so there is no
parent-folder path confusion.

---

## The design

```
100 µm fiber (butted to a 30 µm slit, no relay)
  → f = 78 mm collimator  (scaled Thorlabs AC254-045-B, Ø25.4)
  → Ø16 aperture stop     (f/4.9)
  → 600 l/mm reflective grating  (α = 35.64°, 45° deviation, blaze ~700 nm)
  → f = 78 mm camera      (scaled AC254-045-B, Ø40)
  → plano-cylindrical N-BK7 field flattener  (f = −70 mm, curved in the dispersion plane only)
  → TCD1304 detector (3648 × 8 µm = 29.1 mm), tilted −6°, at 131 mm along the arm
```

- **No fiber relay.** The fiber face is the object at the slit plane (Z = 0). The trace is the
  spectrograph alone, so no relay aberrations enter the resolution figure. A 100 µm fiber
  overfills the 30 µm slit, which keeps the wavelength calibration stable under small alignment
  drift while preserving fine resolution.
- **The key element is the cylindrical field flattener.** The single-doublet camera's limiter
  is *dispersion-plane field curvature* (best focus swings ~2 mm across the band). A cylinder
  curved only in the dispersion plane flattens that focal surface without magnifying the
  slit-height axis — taking the design from ~5.5 nm (uncorrected) to ~1.9 nm.
- Both powered doublets are **one prescription** (AC254-045-B scaled to f = 78), so the
  achromatic correction carries over exactly.

### Scope
This is a **geometric** design validated in BeamFour, and the ~1.9 nm is the single-exposure
delivered figure. Each wavelength is already ~0.5 nm at *its own* best focus, so the optics
themselves aren't the wall — the single flat detector is, together with the small chromatic part
of the focal curve that one cylinder can't fully remove. Closing that to a genuine single-shot
**~1 nm** would take a **field-corrected, multi-element camera** (a Petzval-type group with
enough degrees of freedom for the full field) — a real optical-design job, not a part swap. We
haven't taken that on; ~1.9 nm is the design we're putting forward. The cylindrical flattener is
a **custom** part (like the scaled f = 78 doublet) — both need sourcing or custom fabrication.

---

## Files

| File | Role |
|---|---|
| `med_gen.py` | Sellmeier → `Broadband.MED` (N-LAK22 / N-SF6HT / N-BK7, 400–1000 nm) |
| `frontend_gen.py` | Prescription generator → `Broadband.OPT` (run with **no args** for the released design; cylinder written directly in BeamFour's `Cx` column) |
| `gen_rays.py` | Fiber-cone ray generator → `Broadband.RAY` (auto-caps at BeamFour's 3600-ray limit) |
| `trace_check.py` | Standalone offline ray tracer; reproduces BeamFour to <0.1 µm — for debugging |
| `analyze_spots.py` | BeamFour-traced `Broadband.RAY` → `spots.json` (per-wavelength spot stats) |
| `build_spot_html.py` | `spots.json` → `JASPER_spot_analysis.html` (the report) |
| `Broadband.OPT` | **The file you load into BeamFour** (single prescription; cylinder already in `Cx` form) |
| `Broadband.MED` / `Broadband.RAY` | Glasses and input rays for BeamFour |
| `spots.json` / `JASPER_spot_analysis.html` | Reference traced result and report (~1.9 nm) — to compare against |

---

## Reproduce it

Requires Python 3 + NumPy, and [BeamFour](https://www.stellarsoftware.com/) (free) for the trace.

```bash
pip install -r requirements.txt

# 1. Generate the three BeamFour inputs (Broadband.MED, Broadband.OPT, Broadband.RAY) — one line:
python med_gen.py && python frontend_gen.py && python gen_rays.py spots

# 2. In BeamFour (do this by hand — BeamFour is a GUI):
#    - File → Open  Broadband.OPT     (this replaces the optics table in memory!)
#    - Open Broadband.MED and Broadband.RAY
#    - Run → Layout   (confirm all six wavelengths reach the detector)
#    - Run → InOut    (fills the xfinal/yfinal columns)
#    - SAVE the ray table back to  Broadband.RAY  IN THIS DIRECTORY

# 3. Analyse and report — one line:
python analyze_spots.py && python build_spot_html.py
```

`frontend_gen.py` writes the cylindrical flattener straight into BeamFour's `Cx` column, so
`Broadband.OPT` is the file you load directly — there is no conversion step and only one `.OPT`.

**Sanity check:** `analyze_spots.py` prints the spectrum span. It must read **≈ 28.37 mm**
("CORRECT 131mm design"). If it prints something else, BeamFour re-traced with stale optics —
re-open `Broadband.OPT` and save the trace into this folder. The offline tracer matches BeamFour
exactly, so you can also spot-check geometry with `python trace_check.py`.

### Common gotchas (these cost real time)
- **BeamFour caches the last-loaded OPT.** You must `File → Open Broadband.OPT` *before* InOut,
  or it re-traces with whatever optics were previously in memory.
- **The traced RAY must be saved into this directory** as `Broadband.RAY`. `analyze_spots.py`
  reads it from here; a stale copy elsewhere gives old numbers.
- **BeamFour silently truncates ray files over 3600 rows**, dropping whole wavelengths off the
  end. `gen_rays.py` enforces the cap.
- **The cylinder lives in the `Cx` column** (curved in dispersion, flat along the slit). If the
  traced resolution comes out ~5.5 nm, BeamFour read it in the wrong plane — check that the
  flattener row has its curvature under `Cx` (not `Cy`) and re-trace.

---

## Licensing

Each file is covered by exactly one licence — see this folder's `LICENSE`, and the
repository-root `LICENSE` allocation + `LICENSES/` full texts, which govern:

- **Code** (the `.py` scripts): **MIT** — the tools that generate and analyse the design.
- **Optical-design source** (`.OPT` / `.MED` / `.RAY`, the prescription, this README, and the
  generated `spots.json` / `JASPER_spot_analysis.html`): **CERN-OHL-S-2.0**, strongly-reciprocal
  open hardware, in keeping with JASPER's right-to-repair approach. SPDX: `CERN-OHL-S-2.0`
  (full text: <https://ohwr.org/cern_ohl_s_v2.txt>).

The scripts are permissive; a prescription they *emit* is Covered Source under CERN-OHL-S-2.0.

## Attribution

- Optical ray tracing by **BeamFour** (Stellar Software / M. Lampton et al.).
- JASPER is an open-source VIS-NIR spectrometer platform by **CheckAg**.

Corrections and "have you tried…" welcome via issues.
