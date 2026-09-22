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

| f-number | mean Δλ | worst (500 nm edge) | astig line (spherical) |
|---|---|---|---|
| f/5.3 (wide open) | ~2.24 nm | ~2.74 nm | ~1.2 mm |
| f/8 (toroidal) | ~1.38 nm | ~2.03 nm | ~0.09 mm |
| f/10 + VLS correction | ~0.89 nm | ~1.03 nm | ~0.2 mm |
| aberration-free floor | ~0.6 nm | — | — |

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

## Files

| File | Purpose |
|---|---|
| `vls_raytrace.py` | **The standalone VLS ray tracer — the engine.** Exact sphere/toroid intersection + local VLS grating equation. Imported by every layout/spot script. |
| `lpf_resolution.py` | Light-path-function (analytic) resolution for the design. |
| `lpf_spots.py` | Per-wavelength analytic spot diagrams. |
| `lpf_compare_fig.py` | Analytic resolution / astigmatism comparison figure. |
| `raytrace_two_panels.py` | Ray-traced spots: f/5.3 spherical vs f/8 toroidal. |
| `raytrace_three_panels.py` | Ray-traced spots: the f/5.3 → f/8 → f/10+VLS path to ~1 nm. |
| `raytrace_layout.py` / `_scale.py` / `_light.py` | The three optical-layout renderers. |
| `Noda_Namioka_theory.md` | The aberration theory write-up, with original-paper references. |
| `optiland_experimental/vls_grating.py` | Optiland VLS surface model — validated but **not used** in any delivered result. |
| `figures/` | Pre-generated layouts, spots, comparison. |

## Requirements & run

```bash
pip install -r requirements.txt        # numpy, scipy, matplotlib
```

```bash
python raytrace_layout_light.py        # -> jasper_layout_light.png   (classic traced-ray layout)
python raytrace_layout_scale.py        # -> jasper_layout_toscale.png (to-scale, dimensioned)
python raytrace_layout.py              # -> jasper_raytrace_layout.png (dark)
python raytrace_three_panels.py        # -> jasper_three_panel_spots.png (path to ~1 nm)  [slow: fits VLS]
python raytrace_two_panels.py          # -> jasper_two_panel_spots.png
python lpf_resolution.py               # -> analytic resolution table
```

The delivered pipeline is **pure NumPy + Matplotlib** — Optiland is not required for any result
(only to run `optiland_experimental/vls_grating.py`'s self-test). Every script is deterministic.

## Honesty notes

- **Simulated design**, not a commercial part; **no hardware built**.
- The **f/10 + VLS** ~1 nm figure is **best case** — it assumes the fabricated groove profile matches
  the fitted VLS coefficients. For a real blank the *classical* (uncorrected) numbers are the
  guaranteed floor and the VLS numbers are the plausible target.
- A hard floor sits under everything: the 20 µm slit + pixel + diffraction set the ~0.6 nm
  aberration-free floor. Going lower needs a narrower slit (less light) or more dispersion.

## License

Code (`*.py`) under MIT; optical designs / generated outputs under CERN-OHL-S-2.0, per the parent
repository.
