# JASPER × Shimadzu P0400-01 — buildable spectrometer (Release 3 draft)

A fiber-fed 400–1100 nm spectrometer built around the Shimadzu P0400-01 flat-field grating and catalog
parts. Every number below is from the end-to-end model (`relay_design.py`): real lens surfaces from the
Edmund datasheets → slit → our reconstructed P0400 grating → detector plane. **Modelled, not built.**

## Predicted performance (compact relay, no grating stop)

| | value |
|---|---|
| Band | 400–1100 nm on a 27.66 mm detector (datasheet 27.7 mm) |
| Resolution (FWHM of the traced line) | **1.43 nm mean**; 0.63–0.68 nm near 1000 and 500 nm, up to 2.1 nm at 400 nm |
| Geometric throughput | **~1.05 % of the fiber's light lands on the 200 µm pixels** (0.97–1.08 % across the band) |
| Fiber face → slit | **87 mm** (the 25 mm-lens alternative needs 176 mm) |
| Operating aperture | set by the relay, ≈ f/8–9 at the grating — no grating stop needed |
| Focus tolerance | slit axial position ±1 mm costs ≤ 5 % of the light |

Resolution is quoted as the FWHM of the line profile you would measure on a lamp line (traced, including
the 8 µm pixel and diffraction). The package's quick formula (2.355 × RMS in quadrature) gives 1.20 nm
here; it underestimates wherever defocus dominates, because a defocused line is flat-topped, not Gaussian.
The narrowest window holding 80 % of the line's light (EE80) is 1.27 nm mean (`resolution_metrics.py`).
Throughput is geometric only (étendue, lens aberrations, slit, pixel height); coating, fiber,
grating-efficiency and detector-QE losses come on top.

## Optical train (along the relay axis, from the fiber face)

| Element | Position | Orientation / note |
|---|---|---|
| Fiber face (SMA905) | 0 | 200 µm core, NA 0.22 |
| #49-321 achromat, 14 mm FL, 12.5 mm dia | first vertex at 9.9 mm | **reversed**: the flatter side (R = 56.54 mm) faces the fiber; beam ≈ 6.3 mm dia after it |
| air gap | 5 mm | collimated space |
| #49-329 achromat, 60 mm FL, 12.5 mm dia | first vertex at 22.7 mm | strongly curved side (R = 36.51 mm) faces the collimator |
| 25 µm slit | 87.0 mm (≈ 59.3 mm after the #49-329 rear vertex) | slit length along the grooves; fiber image ≈ 0.86 mm dia |
| P0400-01 grating | slit → pole 96.9 mm, α = 20° | relay axis aimed at the grating pole; mind the blank's orientation mark (rotated 180° → ~10 nm) |
| Detector | datasheet plane: r′ᵢ = 96.9 mm from the pole, normal at βᵢ = −7.51° | ends at r′ = 97.1 mm (400 nm) and 99.5 mm (1100 nm) |

Relay magnification 60/14 = 4.29: NA 0.22 at the fiber → NA 0.051 at the slit.

## Parts (Edmund Optics India prices, incl. taxes, qty 1)

| Part | Stock # | Price |
|---|---|---|
| 0.22 NA VIS/NIR patch cord, 200 µm core, SMA905, 1.5 m (300–2400 nm) | #23-692 | ₹31,382 |
| 12.5 mm dia × 14 mm FL VIS-NIR achromat (N-BAF10 / N-SF5) | #49-321 | ₹9,736 |
| 12.5 mm dia × 60 mm FL VIS-NIR achromat (N-BK7 / N-SF5) | #49-329 | ₹9,736 |
| 25 ±3 µm × 3 mm mounted precision air slit (Ø25 mm) | #58-541 | ₹15,941 |
| **Subtotal, sourced optics** | | **₹66,795** |

Still to source: the P0400-01 grating (Shimadzu quote), the detector board (TCD1304-class, 8 × 200 µm
pixels), mounts for **Ø12.5 mm** lenses (not Ø½" / 12.7 mm — use 12.5 mm-bore mounts or adapters) and the
Ø25 mm slit, an SMA905 fiber adapter, and a light-tight enclosure with baffles.

Alternative (not recommended): #49-352 (30 mm) + #49-361 (125 mm), Ø25 mm — same resolution, ~5 % less
light, twice the length (176 mm), ₹8,172 more.

## Alignment outline

1. Collimate: fiber at the #49-321 back focus (9.9 mm from its flat-side vertex); check with a distant target at ~590 nm.
2. Place the slit ~59.5 mm behind #49-329 and centre it on the fiber image (~0.86 mm across, so lateral
   placement is forgiving); axial ±1 mm is fine.
3. Aim the relay axis at the grating pole, 96.9 mm away at α = 20°, grating orientation mark as specified.
4. Detector on the datasheet plane; fine-tune on a Hg/Ar line lamp.

## Honest notes

- The P0400 groove chirp is our reconstruction from the datasheet (R from the blank sag, chirp solved to
  focus on the published detector plane). It reproduces every published geometric number, but the real
  grooves are Shimadzu's; the measured sample is the ground truth.
- Slit width tolerance ±3 µm ≈ ±0.08 nm at the floor.
- Edmund specifies the lens coating for 400–1000 nm; 1000–1100 nm will see somewhat higher reflection loss.
- The largest remaining loss is vertical: P0400's spherical blank smears each wavelength into a tall
  astigmatic line, so only 30–39 % of the slit light lands on 200 µm-tall pixels. A cylindrical lens in
  front of the detector, or a detector with taller pixels, would recover most of it.
- Why this configuration: `p04_entrance_options.png` — relay + slit at ~f/10 beats using the fiber core as
  the slit, and a 200 µm core gives ~1.8× the light of 100 µm at the same resolution.

Reproduce: `python relay_optics.py` (glass + lens checks) then `python relay_design.py`
(`relay_layouts.png`, `relay_spots_compact.png`, `relay_performance_compact.png`), then `python combined_layout.py` for the full to-scale layout (`jasper_p0400_layout.png`).
