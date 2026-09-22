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

**Analytic resolution (light-path-function / Noda theory)**

| Script | Purpose |
|---|---|
| `lpf_resolution.py` | Resolution for all designs, analytically. **Writes `lpf_all_designs.json`** — run this before `lpf_compare_fig.py`. |
| `lpf_spots.py` | Per-wavelength spot diagrams from the analytic model (best-case). |
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
| 1 | `python vls_raytrace.py` | console only | Reproduces the analytic classical number and 25.3 nm/mm dispersion for P0400 → tracer is sound. (~seconds) |
| 2 | `python lpf_resolution.py` | `lpf_all_designs.json` + table | Simulated rows: plate 34.2 nm/mm, classical floor ~3.7 nm, best-case ~0.6 nm, astig ~1 mm (sph) / ~0.04 mm (tor). |
| 3 | `python lpf_spots.py` | `jasper_simulated_spots.png`, `jasper_p0400_spots.png` | Analytic best-case spots (~0.9–1.05 nm simulated). |
| 4 | `python lpf_compare_fig.py` | `jasper_lpf_comparison.png` | 3 bars: resolution 3.70 / 3.70 / 0.93 nm; astig 1.03 / 0.04 / 0.04 mm. *(needs step 2's JSON)* |
| 5 | `python raytrace_two_panels.py` | `jasper_two_panel_spots.png` | f/5.3 spherical ≈ 2.24 nm mean; f/8 toroid ≈ 1.38 nm. |
| 6 | `python raytrace_three_panels.py` | `jasper_three_panel_spots.png` + `jasper_three_panel_numbers.json` | Panel A 2.24 → B 1.38 → C 0.89 nm mean. **(slow, ~2–3 min: the VLS fit)** |
| 7 | `python raytrace_layout_light.py` | `jasper_layout_light.png` | β 500 nm ≈ −5.7°, 1000 nm ≈ +0.9°; detector tilt ~6°; spectrum ~15 mm. |
| 8 | `python raytrace_layout_scale.py` | `jasper_layout_toscale.png` | L_A 137.36, L_H 126.75, dimensioned. |
| 9 | `python raytrace_layout.py` | `jasper_raytrace_layout.png` | Dark cinematic layout, ~1080 rays. |

**Cross-validation (the "modeled two independent ways that agree" claim):** steps 2 and 5–6 agree
exactly on dispersion and astigmatism; the ray-traced resolution sits at the analytic best-case
floor when stopped down and rises above it at fast f — and that gap **is** the f-number tradeoff.
If those line up, the design is internally consistent.

Steps 1–5 and 7–9 each take seconds to ~1 minute; only step 6 is slow.

---

## Figures to publish (Phase 1)

- `jasper_layout_light.png` — the single-grating concept (real traced rays).
- `jasper_three_panel_spots.png` — the resolution path to ~1 nm (the f-number tradeoff).
- `jasper_lpf_comparison.png` — what each correction buys: VLS fixes resolution, the toroid fixes throughput.

All simulated-only; resolution is presented as the honest f-number tradeoff, not a single number.
