# JASPER simulated VLS grating — scripts & validation order

How to run and validate the simulated concave VLS grating design (230 l/mm, α = 12.34°,
arms 137.357 / 126.75 mm, R = 133 mm, spherical + toroidal), for the 500–1000 nm band.

Everything is **pure NumPy + Matplotlib** — Optiland is not required for any result.

```
pip install -r requirements.txt      # numpy, scipy, matplotlib
```

---

## What each script is

**Engine**

| Script | Purpose |
|---|---|
| `vls_raytrace.py` | The standalone geometric ray tracer (no Optiland). Every ray-traced script imports it. Its own `__main__` runs a self-test on the Shimadzu P0400 grating. |

**Cross-validation (the two models, side by side)**

| Script | Purpose |
|---|---|
| `validate_models.py` | Runs **both** models (analytic LPF + ray tracer) on the same designs and prints a table with a delta column per metric → `jasper_model_validation.png`. This is the machine-checkable form of "modeled two independent ways that agree." **Slow** (runs the VLS fits). |
| `recording_optimize_sim.py` | **Realizability check.** Forward-optimizes a real 2-point holographic recording (real laser lines + geometry, pole-density constrained) and compares the achievable resolution to the free-polynomial best case → `jasper_recording_optimize.png`. Confirms the ~0.9 nm is reachable by a real exposure. |

**Analytic resolution (light-path-function / Noda theory)**

| Script | Purpose |
|---|---|
| `lpf_resolution.py` | Resolution for all designs, analytically. **Writes `lpf_all_designs.json`** — run this before `lpf_compare_fig.py`. |
| `lpf_spots.py` | Per-wavelength spot diagrams from the analytic model — the idealised bound (defocus nulled at every λ). |
| `lpf_compare_fig.py` | The 3-bar comparison: spherical (no VLS) → toroidal (no VLS) → toroidal + VLS. **Reads `lpf_all_designs.json`.** |

**Ray-traced spots (the simulated design)**

| Script | Purpose |
|---|---|
| `raytrace_two_panels.py` | f/5.3 spherical vs f/8 toroidal spots → `jasper_two_panel_spots.png`. |
| `raytrace_three_panels.py` | The f/5.3 → f/8 → f/10+VLS "path to ~1 nm" → `jasper_three_panel_spots.png`. **Slowest** (runs the VLS fit, ~2–3 min). |

**Optical layouts (the simulated design)**

| Script | Purpose |
|---|---|
| `raytrace_layout_light.py` | Light, classic traced-ray layout → `jasper_layout_light.png`. |
| `raytrace_layout_scale.py` | To-scale, dimensioned layout → `jasper_layout_toscale.png`. |
| `raytrace_layout.py` | Dark, cinematic layout → `jasper_raytrace_layout.png`. |

**Kept but not used**

| Script | Purpose |
|---|---|
| `optiland_experimental/vls_grating.py` | Retired Optiland model — kept, **not used** in any result. |

**Theory**

| File | Purpose |
|---|---|
| `Noda_Namioka_theory.md` | The aberration-theory write-up with original-paper references. |

---

## How many rays

Each spot diagram is computed from **2,500 traced rays per wavelength** — across the 6-wavelength
band (500–1000 nm) that is **15,000 rays per panel**. Centroids use 1,500/wavelength; the internal
VLS-fit and best-focus steps use 250/wavelength. Rays are sampled uniformly over the grating
aperture disk (radius EPD/2) with a **fixed random seed**, so every run is deterministic and
reproducible. Layout figures draw far fewer on purpose (15/wavelength for the light layout,
~180/wavelength ≈ 1,080 total for the dark one) — just enough to read the fan.

---

## Execution & validation order

Start by reading and running `vls_raytrace.py` — its self-test is the first proof the tracer is
sound. Then run in this order:

| # | Command | Produces | Check |
|---|---|---|---|
| 1 | `python vls_raytrace.py` | console only | Reproduces the ~25.5 nm/mm dispersion for P0400 exactly → the grating core is sound. (Its classical *resolution* is not comparable to the analytic number: this self-test puts the detector perpendicular to the mid-band ray, while the datasheet detector is tilted. `validate_models.py` makes the like-for-like comparison on the datasheet plane — 8 % agreement.) (~seconds) |
| 2 | `python lpf_resolution.py` | `lpf_all_designs.json` + table | Simulated rows: plate 34.2 nm/mm, classical floor ~3.7 nm, idealised bound ~0.93 nm (not reachable on a flat detector), astig ~1 mm (sph) / ~0.04 mm (tor). |
| 3 | `python lpf_spots.py` | `jasper_simulated_spots.png`, `jasper_p0400_spots.png` | Analytic idealised-bound spots (~0.9–1.05 nm simulated) — labelled as such in the titles. |
| 4 | `python lpf_compare_fig.py` | `jasper_lpf_comparison.png` | 3 bars: resolution 3.70 / 3.70 / 0.93 nm (bar 3 = idealised bound, noted on the figure); astig 1.03 / 0.04 / 0.04 mm. *(needs step 2's JSON)* |
| 5 | `python raytrace_two_panels.py` | `jasper_two_panel_spots.png` | f/5.3 spherical ≈ 2.24 nm mean; f/8 toroid ≈ 1.38 nm. |
| 6 | `python raytrace_three_panels.py` | `jasper_three_panel_spots.png` + `jasper_three_panel_numbers.json` | Panel A 2.24 → B 1.38 → C 0.89 nm mean. Panel C's grooves are the free-polynomial fit (theoretical optimum); step 11 proves a real recording reaches the same ~0.9 nm. **(slow, ~2–3 min: the VLS fit)** |
| 7 | `python raytrace_layout_light.py` | `jasper_layout_light.png` | β 500 nm ≈ −5.7°, 1000 nm ≈ +0.9°; detector tilt ~6°; spectrum ~15 mm. |
| 8 | `python raytrace_layout_scale.py` | `jasper_layout_toscale.png` | L_A 137.36, L_H 126.75, dimensioned. |
| 9 | `python raytrace_layout.py` | `jasper_raytrace_layout.png` | Dark cinematic layout, ~1080 rays. |
| 10 | `python validate_models.py` | `jasper_model_validation.png` + `model_validation.json` | LPF vs ray trace, delta per metric (see below). **(slow: VLS fits)** |
| 11 | `python recording_optimize_sim.py` | `jasper_recording_optimize.png`, `jasper_recording_spots.png` + `recording_optimize_sim.json` | Best real recording ≈ 0.90 nm mean (matches the 0.88 nm free polynomial → realizable); per-λ spots for the physical grooves. **(slow: VLS fit + laser scan)** |
| 12 | `python resolution_metrics.py` | `jasper_resolution_metrics.png` + `resolution_metrics.json` | Quick formula vs measured line widths (finite slit traced; FWHM and EE80 = window holding 80 % of the light). Panel C f/10: formula 0.89 → FWHM 0.76 / EE80 0.74 nm (worst 0.93 / 0.94). 488 nm recording: 0.90 → 0.74 / 0.75. Panel A f/5.3: 2.24 → FWHM 0.94 but EE80 2.01 (coma tail). P0400 build: 1.20 → FWHM 1.43 / EE80 1.27 (flat-topped). **(slow, ~7 min; needs `validate_p0400.json` and `recording_optimize_sim.json` in the working directory)** |

**Two groove models, kept separate.** `vls_raytrace.py` represents the correction as a *polynomial*
ν(y) = 230 + c₁y + c₂y² + c₃y³ (the c's fitted by `fit_vls`). `recording_optimize_sim.py` represents
it as the *physical* holographic function ν(P) = (1/λ₀)(û_C − û_D)·ĝ from two recording sources —
applied **directly** in its own tracer, **not** converted to c₁,c₂,c₃ and passed to `vls_raytrace`.
The polynomial path is the reference (via `fit_vls`/`evaluate`); the holographic path is the
realizable design. The script prints both and confirms they land on the same ~0.9 nm.

**Cross-validation (the "modelled two independent ways" claim):** run `python validate_models.py`.
Both models are evaluated on the **identical fixed flat detector plane** (normal along β_ref, same
distance, no refocusing), and every metric is like-for-like. Results (f/5.3, mean over the band):

| | Sim spherical | Sim toroidal | P0400-01 (datasheet geometry) |
|---|---|---|---|
| plate scale | 0.1 % | 0.1 % | 0.0 % |
| astigmatic line | 0.1 % | 0.04 vs 0.10 mm (both collapsed) | 0.1 % |
| classical resolution (constant grooves) | **1.3 %** | 13 % | **8 %** |
| idealised bound (defocus nulled at every λ) | 0.93 vs 0.75 nm | 0.93 vs 1.12 nm | **0.613 vs 0.611 nm** |
| one real chirp (identical c₁ in both) | 2.23 vs 2.55 nm | 1.87 vs 2.61 nm | **2.36 vs 2.72 nm** |

- **Classical** is the cleanest check: same plane, constant grooves, no free knobs.
- **Idealised bound**: both models assume defocus is nulled separately at every wavelength (the
  analytic model by setting F₂₀ = 0, the ray tracer by fitting a separate chirp per wavelength). It is a
  *bound*, not a design: no real grating reaches it on a flat detector, because a single groove chirp —
  a VLS polynomial *or* a two-point holographic recording — adds a defocus correction ∝ λ × constant and
  so nulls defocus at only ~2 wavelengths. (The ray tracer's per-wavelength fit also corrects coma,
  which is why it can sit below the analytic value.)
- **One real chirp**: the physically buildable case. One linear chirp c₁ is fitted once in the ray tracer
  and fed unchanged to the analytic model. Both show the same W-shaped curve — floor at the two defocus
  nulls, rising between and beyond. Letting the ray tracer also use c₂ (coma) lowers P0400 to 2.02 nm
  mean, 0.68 nm at the null.
- The jump from ~0.6 nm (idealised) to ~2–2.7 nm (one chirp) is **not** the two models disagreeing —
  both show it. It is the flat-field compromise itself.
- Gaps for the toroid are larger (13–40 %) because the analytic formulas keep only the meridional terms;
  the off-meridian term that bends the spectral lines depends on a near-cancellation (1/r − cos α/ρ) that is
  very sensitive to the sagittal radius ρ. The script's mid-band sign check shows it: the single-wavelength
  c₁ matches the analytic value for both spherical gratings (0.0383 vs 0.0385; 0.1333 vs 0.1338) but not
  for the toroid (0.024 vs 0.039).
- For the **simulated** design the detector position is a free design choice; this check pins it to the
  analytic model's plane only to keep the comparison like-for-like. The design figures
  (`raytrace_three_panels.py`) optimise the detector distance and stop down to f/10, which is where the
  ~0.9 nm comes from — always a ray-traced number, never the idealised bound (measured on the traced lines
  it is ~0.75 nm; see step 12).

*Correction (how resolution is measured):* earlier numbers used a quick formula, hypot(slit, 2.355 × RMS,
pixel, diffraction) × plate scale, which assumes bell-shaped lines. `resolution_metrics.py` now measures the
traced lines directly (finite slit traced, 8 µm pixel, diffraction) and reports FWHM and EE80. The formula
errs in both directions: it reads ~16 % narrow for the P0400 build (flat-topped, defocus-limited lines:
1.20 vs 1.43 nm FWHM) and ~15 % wide for the Release 1 f/10 VLS design (compact lines: 0.89 vs 0.76 nm FWHM,
0.74 nm EE80). For lines with a coma tail (uncorrected, f/5.3) FWHM alone misleads — 0.94 nm FWHM while
80 % of the light spans 2.01 nm — so both widths are reported, and EE80 is used to compare f-numbers.

*Correction (plate scale, after Shimadzu's brochure eq. 24):* on a flat detector tilted to the
β = βᵢ ray, the dispersion along the detector is D = cos β · cos(βᵢ − β) / (N m r′), with r′ the
distance along the λ ray (r′ = r′ᵢ / cos(βᵢ − β), eq. 23). Earlier versions of the analytic plate scale
used cos β / r′ (missing one cos(βᵢ − β)) and the P0400 relay/entrance scripts used
cos β / (r′ cos(βᵢ − β)) (inverted). Both overstated nm/mm toward the ends of the detector — by up to
5.3 % at 1100 nm for P0400 — so those resolutions were slightly pessimistic. Every script now uses
eq. 24; the ray-traced numbers (plate from traced centroids) were always correct and did not change.
Effect: LPF-vs-trace plate agreement 0.8 % → 0.0 % for P0400; relay FWHM 1.44 → 1.42 nm mean;
validate_p0400 (2.02 / 1.19 nm) unchanged.

*Correction (detector plane, after the datasheet table):* the P0400 plane was first fitted through the
rounded r′₁ = 97.1 / r′₂ = 99.5 mm (96.96 mm, −7.36°). The datasheet also prints the plane itself,
r′ᵢ = 96.9 mm at βᵢ = −7.51°; with that βᵢ, all four printed values agree for r′ᵢ = 96.92–96.95 mm. Every
script now uses βᵢ = −7.51°, r′ᵢ = 96.93 mm. Effect: C₁ +0.2539 → +0.2530 (f/5.3), detector 27.65 → 27.66 mm,
resolution 2.02 / 1.19 nm unchanged, relay FWHM 1.42 → 1.43 nm, slit 0.2 mm closer to the relay.

*Correction:* an earlier version of this check reported a ~60 % classical gap for P0400 and attributed
it to higher-order aberration. Most of it was ours: the ray tracer and the analytic model were using
different detector planes. On the same plane the gap is 8 %.

Steps 1–5 and 7–9 each take seconds to ~1 minute; step 6 and `validate_models.py` are slow (VLS fits).

---

## Release 3 — Shimadzu P0400-01 validation and a buildable spectrometer

Validates the pipeline on a real commercial grating, then designs a spectrometer around it from
catalog parts. Run in this order (each writes its figure + JSON; later steps read earlier JSON):

| # | Script | Produces | Check / result |
|---|---|---|---|
| R1 | `python validate_p0400.py` | `validate_p0400.png`, `p0400_groove_profile.png`, `validate_p0400.json` | Datasheet reproduced: β −10.49°/+5.62°, detector 27.66 mm (27.7); the detector plane is taken from the published r′ᵢ / βᵢ (an input, not a check). R = 99.61 mm from the blank sag (free fit 100.00). Solved chirp c₁ ≈ +0.253; 400 l/mm at the pole, peak 400.9 at +7.4 mm. Resolution 2.02 nm (f/5.3) / 1.19 nm (f/10). |
| R2 | `python p0400_verify.py` | `p0400_verify.json` | Pole vs edge hypothesis test (pole pinned to ±0.16 l/mm); blank rotated 180° → 9.95 nm; idealised bound vs single chirp in both models. |
| R3 | `python p0400_recording.py` | `p0400_recording.png`, `.json` | Two-point holographic recordings for 7 laser lines all reproduce the chirp (2.01–2.02 nm); λ₀ not identifiable from optical data; two source-side families differ ~1.3 l/mm at the blank edges. |
| R4 | `python p0400_entrance_options.py` | `p04_entrance_options.png`, `p0400_entrance_options.json` | Fiber → spectrometer entrance options (fiber-as-slit, relay + slit, slit butted to fiber; 50–400 µm cores; f/5.3 and f/10): resolution vs light reaching the 200 µm pixels. Relay + slit at f/10 wins; 200 µm core ≈ 1.8× the light of 100 µm at the same resolution. |
| R5 | `python relay_optics.py` | console | Glass check (Sellmeier vs catalog nd/vd, exact) and lens check (Edmund prescriptions reproduce EFL/BFL within 0.04 mm). |
| R6 | `python relay_design.py` | `relay_layouts.png`, `relay_spots_compact.png`, `relay_performance_compact.png`, `relay_design_<name>.json` | End-to-end real-surface trace, fiber #23-692 → relay → 25 µm slit #58-541 → P0400 → pixels, for two catalog relays. Recommended compact relay #49-321 (14 mm) + #49-329 (60 mm): fiber → slit 87 mm (vs 176 mm for the 25 mm-lens pair), 1.05 % of the fiber's light on the pixels, resolution 1.43 nm mean as measured FWHM of the traced line (0.63–2.12 nm); the quick Gaussian formula gives 1.20 nm and underestimates where defocus dominates. |
| R7 | `python combined_layout.py`, `python combined_layout_dark.py` | `jasper_p0400_layout.png` (light), `jasper_p0400_layout_dark.png` (dark) | Complete optical layout to scale with traced rays: fiber → compact relay → slit → P0400-01 → detector, real surfaces throughout; reports the slit-mount-to-detector clearance (~11 mm to the active area). Needs R6's `relay_design_compact.json`. |

Supporting modules: `relay_optics.py` (glass dispersion, catalog prescriptions, real-surface 3-D ray tracer).
R1 is slow (~2 min, chirp fits); R3–R5 take ~1–2 min each, R6 ~4 min (both relays). See `JASPER_P0400_build.md` for the parts list.

## Figures to publish (Phase 1)

- `jasper_layout_light.png` — the single-grating concept (real traced rays).
- `jasper_three_panel_spots.png` — the resolution path to ~1 nm (the f-number tradeoff).
- `jasper_lpf_comparison.png` — what each correction buys: VLS fixes resolution, the toroid fixes throughput.
- `jasper_resolution_metrics.png` — quick formula vs measured FWHM and EE80 for every headline design.

All simulated-only; resolution is presented as the honest f-number tradeoff, not a single number.
