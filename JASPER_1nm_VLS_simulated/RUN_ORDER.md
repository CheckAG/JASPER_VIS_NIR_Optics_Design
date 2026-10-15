# Run order

Every script runs from this folder and writes its PNG/JSON outputs here (the shipped copies are in `figures/`).
Setup once:

```
pip install -r requirements.txt
```

Run everything in order with `python run_all.py`, or resume with `python run_all.py --from 12`.
Times are from a clean unzip on a cloud Linux machine; expect 1–2× longer on a laptop.

## Release 1 — simulated custom VLS grating (230 l/mm, 500–1000 nm)

| # | Command | Needs | Writes | Time |
|---|---|---|---|---|
| 1 | `python vls_raytrace.py` | — | console self-test | 15 s |
| 2 | `python lpf_resolution.py` | — | `lpf_all_designs.json` | <1 s |
| 3 | `python lpf_spots.py` | — | `jasper_simulated_spots.png`, `jasper_p0400_spots.png` | 2 s |
| 4 | `python lpf_compare_fig.py` | 2 | `jasper_lpf_comparison.png` | 1 s |
| 5 | `python raytrace_two_panels.py` | — | `jasper_two_panel_spots.png` | 3 s |
| 6 | `python raytrace_three_panels.py` | — | `jasper_three_panel_spots.png`, `jasper_three_panel_numbers.json` | ~2 min |
| 7 | `python raytrace_layout_light.py` | — | `jasper_layout_light.png` | 2 s |
| 8 | `python raytrace_layout_scale.py` | — | `jasper_layout_toscale.png` | 2 s |
| 9 | `python raytrace_layout.py` | — | `jasper_raytrace_layout.png` | 2 s |

## Realizability and model cross-check

| # | Command | Needs | Writes | Time |
|---|---|---|---|---|
| 10 | `python recording_optimize_sim.py` | — | `jasper_recording_optimize.png`, `jasper_recording_spots.png`, `recording_optimize_sim.json` | 30 s |
| 11 | `python validate_models.py` | — | `jasper_model_validation.png`, `model_validation.json` | 35 s |

## Release 3 — Shimadzu P0400-01 practical build

| # | Command | Needs | Writes | Time |
|---|---|---|---|---|
| 12 | `python validate_p0400.py` | — | `validate_p0400.png`, `p0400_groove_profile.png`, **`validate_p0400.json`** | 5 s |
| 13 | `python p0400_verify.py` | — | `p0400_verify.json` | 2 s |
| 14 | `python p0400_recording.py` | — | `p0400_recording.png`, `p0400_recording.json` | 40 s |
| 15 | `python p0400_entrance_options.py` | 12 | `p04_entrance_options.png`, `p0400_entrance_options.json` | 1 s |
| 16 | `python relay_optics.py` | — | console: glass and catalog-lens check | <1 s |
| 17 | `python relay_design.py` | 12 | `relay_spots_compact.png`, `relay_layouts.png`, `relay_performance_*.png`, **`relay_design_compact.json`** | 15 s |
| 18 | `python fiber_na_compare.py` | 12 | console table | 2 s |
| 19 | `python combined_layout.py` | 12, 17 | `jasper_p0400_layout.png` | 1 s |
| 20 | `python combined_layout_dark.py` | 12, 17 | `jasper_p0400_layout_dark.png` | 3 s |

## Resolution re-check (all designs)

| # | Command | Needs | Writes | Time |
|---|---|---|---|---|
| 21 | `python resolution_metrics.py` | 10, 12 | `jasper_resolution_metrics.png`, `resolution_metrics.json` | ~4.5 min |

## Dependencies in one line

`validate_p0400.json` (step 12) feeds 15, 17, 18, 21. `relay_design_compact.json` (17) feeds 19, 20.
`lpf_all_designs.json` (2) feeds 4. `recording_optimize_sim.json` (10) feeds 21; if it is missing, 21 falls
back to the copy in `figures/`. Every other script stands alone. Total ≈ 10 min.

## Headline numbers you should see

| Step | Line to check | Expected |
|---|---|---|
| 6 | `Panel C ... mean` | 0.89 nm (quick formula) |
| 12 | `[f/5.3] ... mean` / `[f/10] ... mean` | 2.02 nm / 1.19 nm |
| 16 | `#49-321: EFL` | 14.01 (spec 14.00) |
| 17 | `compact ... dlam ... (direct FWHM)` | 1.43 nm |
| 21 | `Release 1 panel C ... FWHM` | 0.76 nm mean |
