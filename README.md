# JASPER_VIS_NIR_Optics_Design
VIS-NIR Spectrometer design using Beam Four and design files generated using AI generated python scripts  


# JASPER — VIS-NIR Optics Design

Open, reproducible optical designs for **JASPER**, an open-source, right-to-repair VIS-NIR
spectrometer platform by [CheckAg](https://hackaday.io/project/202421-jasper-vis-nir-spectrometer).

This repository holds the **optical design** for JASPER's spectrometer variants — the ray-trace
models, the Python that generates and analyses them, and everything needed to reproduce each
result from scratch. The designs are done in **[BeamFour](https://www.stellarsoftware.com/)**,
a free geometric ray tracer, so anyone can check them without expensive optics software.

---

## Flagship result

**500–1000 nm, ~1.9 nm mean / 2.9 nm worst, single exposure — validated in BeamFour.**

Built from stock 1-inch optics plus one custom cylindrical field flattener that corrects the
camera's dispersion-plane field curvature. The offline Python tracer reproduces BeamFour to
better than 0.1 µm per wavelength, so the number is checked in two independent tools.

→ Files and full write-up: [`500_1000nm/`](https://github.com/CheckAG/JASPER_VIS_NIR_Optics_Design/tree/main/500_1000nm)

---

## Repository layout

Each **design line** lives in its own folder, self-contained, with its own README and a
reference trace you can compare against:

```
JASPER_VIS_NIR_Optics_Design/
├── README.md            ← you are here
├── LICENSE              ← licence allocation (which file → which licence)
├── LICENSES/            ← full texts: MIT.txt, CERN-OHL-S-2.0.txt
└── 500_1000nm/          ← 500–1000 nm design (~1.9 nm, BeamFour-validated)
    ├── README.md            reproduce steps + gotchas
    ├── *.py                 generators, tracer, analysis, report
    ├── Broadband.OPT      the file you load into BeamFour (cylinder as Cx)
    ├── Broadband.MED/.RAY    glasses + input rays
    └── spots.json / *.html   reference traced result + report
```

| Design line | Range | Detector | Resolution | Status |
|---|---|---|---|---|
| [`500_1000nm/`](https://github.com/CheckAG/JASPER_VIS_NIR_Optics_Design/tree/main/500_1000nm) | 500–1000 nm | TCD1304 (29.1 mm) | ~1.9 nm mean / 2.9 nm worst | BeamFour-validated |

(More lines — e.g. extended-range and alternate-detector variants — will be added as folders.)

---

## How the designs work

Every design follows the same reproducible flow, all inside its folder:

1. **Python generates the BeamFour inputs** — the prescription (`.OPT`), glasses (`.MED`), and
   fiber-launch rays (`.RAY`).
2. **You trace it in BeamFour** (Run → InOut) and save the traced rays.
3. **Python turns the trace into a spot-analysis report** (per-wavelength resolution).

A standalone offline tracer (`trace_check.py`) reproduces BeamFour's geometry independently, so
the design loop is fast and every headline number is cross-checked. See each folder's README for
the exact commands and the BeamFour steps.

**Requirements:** Python 3 + NumPy, and BeamFour (free) for the authoritative trace.

---

## Design philosophy

- **Reproducible and free.** Stock parts where possible, free software, and files anyone can
  re-trace — no black boxes.
- **Honest numbers.** Resolution figures are the delivered (slit ⊕ optics ⊕ pixel) values from
  an actual trace, with custom parts and limitations stated plainly.
- **Right to repair.** Part of JASPER's open-core approach: open hardware and optics, so the
  instrument can be understood, built, and fixed by its owner.

---

## Licensing

This repo mixes software and open-hardware source; each file is covered by exactly one
licence. The repository-root [`LICENSE`](LICENSE) is the authoritative allocation, with full
texts in [`LICENSES/`](LICENSES/):

- **Code** (`*.py`): **MIT** — the tools that generate and analyse the design.
- **Optical-design source** (`*.OPT` / `*.MED` / `*.RAY`, the prescriptions, the design docs,
  and generated outputs like `spots.json` / the report): **CERN-OHL-S-2.0**, strongly-reciprocal
  open hardware. SPDX: `CERN-OHL-S-2.0` · full text: <https://ohwr.org/cern_ohl_s_v2.txt>

The scripts are permissive, but a prescription they *emit* is Covered Source under
CERN-OHL-S-2.0 — modify a script freely; a new design it produces stays open.

## About

JASPER is a configurable VIS-NIR grating spectrometer platform targeting food quality, soil
health, dairy, and environmental sensing. Build logs and background are on
[Hackaday](https://hackaday.io/project/202421-jasper-vis-nir-spectrometer).

Optical ray tracing by BeamFour (Stellar Software / M. Lampton et al.).

Corrections, questions, and "have you tried…" are welcome via
[issues](https://github.com/CheckAG/JASPER_VIS_NIR_Optics_Design/issues).


