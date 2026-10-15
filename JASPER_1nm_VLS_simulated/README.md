# JASPER — simulated ~1 nm concave VLS grating (500–1000 nm)

Part of **[JASPER_VIS_NIR_Optics_Design](https://github.com/CheckAG/JASPER_VIS_NIR_Optics_Design)** —
CheckAg's open-source, right-to-repair VIS–NIR spectrometer platform.

This package is the optical core of JASPER: a **single concave variable-line-spacing (VLS) grating**
that both disperses *and* focuses 500–1000 nm onto a cheap linear CCD — no collimating or camera
mirrors. It is a **simulated / theoretical design**, modeled two independent ways that agree, with
real geometric ray-traced layouts.

> **Status: modeled, not built.** No hardware yet, and this is a theoretical design, not a specific
> commercial part. Every number here is from the model; "best case" numbers additionally assume the
> groove correction matches the design (see *Honesty notes*).

> **New — a buildable version (Release 3).** The theoretical grating below needs a custom holographic
> master. The package now also contains a practical spectrometer built around a catalog grating, the
> Shimadzu P0400-01, fed by a catalog fiber and relay — see
> [*Practical build: Shimadzu P0400-01*](#practical-build-shimadzu-p0400-01-release-3).

## The design

A variable-line-spacing concave holographic grating in a finite-arm, Shimadzu-style self-imaging
mount: the slit sits above the grating normal at the entrance arm, the spectrum images below the
normal at the exit arm on a flat, tilted detector.

| Parameter | Value |
|---|---|
| Grooves | 230 l/mm at the pole, m = 1 |
| Substrate | concave, R = 133 mm (spherical **or** toroidal) |
| Incidence α | 12.34° |
| Entrance / exit arm | rA = 137.36 mm / rB = 126.75 mm |
| Plate scale | 34.2 nm/mm |
| Spectrum on array | ~14.6 mm on the 29.18 mm TCD1304 |
| Detector | TCD1304, 3648 × 8 µm, flat, tilted ~6° |
| Slit | 20 µm × 1 mm, fed by an SMA905 fiber + relay |

## Resolution is an f-number tradeoff — not a single number

The honest headline is a curve, not a figure. Faster f = more light but more higher-order
aberration; slower f = sharper but dimmer. Ray-traced across 500–1000 nm:

| f-number | EE80, mean (worst) | FWHM, mean (worst) | quick formula, mean (worst) | astig line |
|---|---|---|---|---|
| f/5.3 (wide open, spherical) | **2.01 nm** (2.65) | 0.94 nm (1.02) | 2.24 nm (2.74) | ~1.2 mm |
| f/8 (toroidal) | **1.32 nm** (2.11) | 0.99 nm (1.54) | 1.38 nm (2.03) | ~0.09 mm |
| f/10 + VLS correction | **0.74 nm** (0.94) | 0.76 nm (0.93) | 0.89 nm (1.03) | ~0.2 mm |
| aberration-free floor | ~0.6 nm | — | — | — |

*How resolution is measured* (`resolution_metrics.py`, `figures/jasper_resolution_metrics.png`): the finite
slit is traced, the 8 µm pixel and diffraction included, and two widths are read off each traced line —
**FWHM** (what a lamp line shows) and **EE80** (the narrowest window holding 80 % of the line's light). The
first release quoted a quick formula (2.355 × RMS, in quadrature), which assumes bell-shaped lines. At f/10
the line is compact and the formula was conservative: ~0.9 nm quoted, ~0.75 nm measured. Without the VLS
correction the lines carry a long coma tail — a sharp core (FWHM ~1 nm) with most of the light spread over
~2 nm — so **EE80 is the column to compare across f-numbers**; FWHM alone would make f/5.3 look as sharp
as f/10.

**The binding constraint is the 500 nm blue edge, and it is aberration-limited** (at f/8 the 20 µm
slit contributes a flat ~0.69 nm everywhere, but the aberration term at 500 nm is ~1.9 nm). Reaching
a true ~1 nm across the whole band takes three things *together*, none sufficient alone: stop down to
~f/10, keep the toroid (throughput only), and rely on a real VLS groove correction to null the field
aberration at the edges. `raytrace_three_panels.py` traces that A→B→C path.

**Throughput / substrate.** Each wavelength images as a vertical astigmatic line. A **spherical**
blank leaves a ~1 mm line (~19 % of the light on the 200 µm pixel); a **toroidal** blank collapses it
to ~0.04 mm (~100 %). Same resolution, ~10× the light — the toroid is the throughput upgrade, not a
resolution one.

## Modeled two independent ways that agree

1. **Light-path-function (Noda–Namioka) aberration theory** — `lpf_resolution.py`. The standard
   analytic method for concave gratings: the grating equation sets dispersion; the F₂₀ / F₃₀ / F₀₂
   coefficients set the spot sizes from the mount geometry alone. See `Noda_Namioka_theory.md`.
2. **A standalone geometric ray tracer** — `vls_raytrace.py`, pure NumPy. Traces thousands of real
   rays through the exact concave surface with the **local VLS grating equation**.

They agree exactly on dispersion and astigmatism. On resolution, the ray trace additionally captures
the higher-order aberration the 3rd-order analytic theory truncates — which is why the honest number
is the f-number curve above, not the analytic floor alone.

**Run `validate_models.py` for the automated side-by-side** (`figures/jasper_model_validation.png`).
On an identical fixed detector plane the two models agree to ≤0.1 % on dispersion, ~0.1 % on
astigmatism, and 1.3 % on classical resolution for this design (8 % for the real 400 l/mm Shimadzu
P0400-01). They also agree on the flat-field compromise itself: both put the idealised bound (defocus
nulled at every wavelength) near the ~0.6 nm floor, and both put a single real groove chirp at ~2–2.7 nm
mean at f/5.3 — floor at two wavelengths, rising between. The quoted resolution is always the ray-traced
one.

### Why not Optiland or Beam4?

We started in **Optiland** but found its sequential trace **unreliable for this off-axis reflective
mount** — it could match a known dispersion *or* focus, but not both at once, so its spot sizes were
fitting artifacts. **Beam4** renders beautiful ray fans but its grating model is constant-line only —
it cannot represent variable line spacing, which is the whole point of this grating. So the delivered
pipeline is the analytic model + the standalone tracer, both pure NumPy. The Optiland VLS surface
model that *did* validate (constant grooves reproduce the grating equation to ~1e-16) is kept, unused,
in `optiland_experimental/` for anyone pursuing the ray-trace path.

## Optical layouts (real traced rays, three styles)

All three are genuine ray traces from `vls_raytrace.py` — slit → exact concave surface → local VLS
diffraction → tilted detector — not schematics.

| Script | Figure | Use |
|---|---|---|
| `raytrace_layout.py` | `figures/jasper_raytrace_layout.png` | dark, cinematic — web / Hackaday |
| `raytrace_layout_scale.py` | `figures/jasper_layout_toscale.png` | to-scale, dimensioned — engineering / design record |
| `raytrace_layout_light.py` | `figures/jasper_layout_light.png` | light, classic traced-ray — print / paper |

## Practical build: Shimadzu P0400-01 (Release 3)

The theoretical grating above needs a custom holographic master: we showed a two-point exposure can write
it (`recording_optimize_sim.py`), but nobody sells it. This section designs a spectrometer from parts you
can order: the Shimadzu **P0400-01** flat-field grating (400 l/mm, 400–1100 nm), an Edmund 200 µm / NA 0.22
fiber, a two-lens catalog relay and a 25 µm slit.

**The one unknown is the grooves.** Shimadzu publishes the mount but not the groove pattern, so we
reconstruct it from the datasheet — pole density from the band-edge angles, detector plane from the published
mounting parameters (r′ᵢ, βᵢ), radius of curvature from the blank's sag — and solve the chirp
ν(y) = 400 + C₁y + C₂y² + C₃y³ that focuses the band onto that plane: **C₁ = +0.253 l/mm², C₂ = −0.0172 l/mm³,
C₃ ≈ 0**. The reconstruction reproduces the datasheet's band-edge angles exactly and its detector length to 0.04 mm
(the detector plane itself is taken from the datasheet, so it is an input, not a check). Those confirm the pole
density only: the chirp focuses each wavelength but barely moves it (detector length changes < 0.03 mm with or
without C₁–C₃), so C₁–C₃ still need a measurement. Full derivation: [`P0400_VLS_derivation.md`](P0400_VLS_derivation.md).

| Predicted (modelled, not built) | |
|---|---|
| Resolution, FWHM of the traced line | 1.43 nm mean over 400–1100 nm; 0.63–0.68 nm at 1000 and 500 nm; 2.1 nm at 400 nm |
| Light on the 8 × 200 µm pixels | ~1 % of the fiber's light (geometric) |
| Fiber → slit | 87 mm (relay #49-321 + #49-329) |
| Sourced optics | ₹66,795 (fiber, two lenses, slit) + grating + detector |

Resolution here is the **FWHM of the traced line profile** — what a lamp line would show (EE80, the window
holding 80 % of the light: 1.27 nm mean). The P0400's lines are flat-topped because defocus dominates, so the
quick formula (2.355 × RMS in quadrature) reads narrow here: 1.20 nm. For the theoretical design above it
erred the other way. See `resolution_metrics.py`.

| Script (run in order) | Produces |
|---|---|
| `validate_p0400.py` | Datasheet reproduction, R from the blank, solved chirp → `validate_p0400.png`, `p0400_groove_profile.png` |
| `p0400_verify.py` | Pole-vs-edge test, orientation test, idealised vs real chirp in both models |
| `p0400_recording.py` | Two-point holographic recordings that write the chirp, 7 laser lines → `p0400_recording.png` |
| `p0400_entrance_options.py` | Fiber-as-slit vs relay + slit, core sizes, f-numbers → `p04_entrance_options.png` |
| `fiber_na_compare.py` | Does a 0.10 NA fiber help? (no — ~4× less light, same resolution) |
| `relay_optics.py` | Glass data + Edmund prescriptions (verified against catalog nd/vd and EFL/BFL) + 3-D real-surface tracer |
| `relay_design.py` | End-to-end fiber → relay → slit → grating → pixels, two relays → `relay_spots_compact.png`, `relay_layouts.png` |
| `combined_layout.py` / `combined_layout_dark.py` | The complete layout, to scale → `jasper_p0400_layout.png`, `jasper_p0400_layout_dark.png` |

Parts, positions, orientations and alignment: [`JASPER_P0400_build.md`](JASPER_P0400_build.md).

## Files

| File | Purpose |
|---|---|
| `vls_raytrace.py` | **The standalone VLS ray tracer — the engine.** Exact sphere/toroid intersection + local VLS grating equation. Imported by every layout/spot script. |
| `lpf_resolution.py` | Light-path-function (analytic) resolution for the design. |
| `validate_models.py` | **Runs both models side by side** and prints a delta-per-metric table → `figures/jasper_model_validation.png`. The automated form of the "modeled two ways" claim. |
| `recording_optimize_sim.py` | **Forward-optimizes a real 2-point holographic recording** (real laser lines + geometry, pole-density constrained) and checks it against the free-polynomial best case → `figures/jasper_recording_optimize.png`. Shows the ~0.9 nm is physically realizable, not a free-polynomial artifact. |
| `resolution_metrics.py` | **Resolution re-check**: quick formula vs measured FWHM and EE80 (finite slit traced) for every headline design → `figures/jasper_resolution_metrics.png`. |
| `lpf_spots.py` | Per-wavelength analytic spot diagrams. |
| `lpf_compare_fig.py` | Analytic resolution / astigmatism comparison figure. |
| `raytrace_two_panels.py` | Ray-traced spots: f/5.3 spherical vs f/8 toroidal. |
| `raytrace_three_panels.py` | Ray-traced spots: the f/5.3 → f/8 → f/10+VLS path to ~1 nm. |
| `raytrace_layout.py` / `_scale.py` / `_light.py` | The three optical-layout renderers. |
| `Noda_Namioka_theory.md` | The aberration theory write-up, with original-paper references. |
| `VALIDATION.md` | Every check, step by step, with the numbers you should reproduce and the corrections log. |
| `JASPER_VLS_coefficients_derivation.md` | How C₁–C₃ follow from a two-point holographic recording (the simulated grating). |
| `P0400_VLS_derivation.md` / `JASPER_P0400_build.md` | Release 3: C₁–C₃ reconstructed for the P0400-01 / parts, positions and alignment for the build. |
| `optiland_experimental/vls_grating.py` | Optiland VLS surface model — validated but **not used** in any delivered result. |
| `run_all.py` / `RUN_ORDER.md` | Runs every script in dependency order / the documented sequence with expected numbers. |
| `figures/` | Pre-generated layouts, spots, comparison. |

## Requirements & run

```bash
pip install -r requirements.txt        # numpy, scipy, matplotlib
```

```bash
python run_all.py                      # all 21 scripts in dependency order (~10 min)
python run_all.py --from 12            # resume from a step;  --only 12 17  runs just those
```

The full sequence — what each step needs, writes, how long it takes, and the headline number to check —
is in [`RUN_ORDER.md`](RUN_ORDER.md). Outputs are written to this folder; the shipped copies are in `figures/`.

The delivered pipeline is **pure NumPy + Matplotlib** — Optiland is not required for any result
(only to run `optiland_experimental/vls_grating.py`'s self-test). Every script is deterministic.

## Honesty notes

- **Simulated design**, not a commercial part; **no hardware built**.
- The **f/10 + VLS** ~1 nm figure is **best case** — it assumes the fabricated groove profile matches
  the fitted VLS coefficients. For a real blank the *classical* (uncorrected) numbers are the
  guaranteed floor and the VLS numbers are the plausible target.
- **The best case is realizability-checked.** `recording_optimize_sim.py` searches real holographic
  recording-laser lines and geometry (a two-point exposure, not free coefficients) and finds one that
  reaches **0.90 nm mean** — matching the 0.88 nm free-polynomial fit (quick formula; measured on the traced
  lines: 0.74 nm FWHM, 0.75 nm EE80, worst 0.89 nm). So the ~0.9 nm target is
  achievable by a *real* recording, not just an unconstrained-polynomial artifact. A simple two-point
  exposure buys only low-order correction; if a design ever needs more, that means aspheric-wavefront
  or computer-generated-hologram recording — more capable, more expensive.
- A hard floor sits under everything: the 20 µm slit + pixel + diffraction set the ~0.6 nm
  aberration-free floor. Going lower needs a narrower slit (less light) or more dispersion.

## License

Code (`*.py`) under MIT; optical designs / generated outputs under CERN-OHL-S-2.0, per the parent
repository.
