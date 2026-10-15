# Deriving the Shimadzu P0400-01 groove chirp (C₁, C₂, C₃) from its datasheet

Shimadzu publishes everything about the P0400-01 **except its grooves**. It is an aberration-corrected
(holographic, variable-groove-density) flat-field grating, and the groove pattern is proprietary. This note
shows how we reconstruct it from the published numbers, what each coefficient does, how well the result
reproduces the datasheet, and how to measure it on a real sample.

The groove density across the blank is written as

$$\nu(y) = \nu_0 + C_1 y + C_2 y^2 + C_3 y^3 \qquad [\text{lines/mm}]$$

with *y* measured from the grating pole (centre), positive toward the entrance-slit side.

Scripts: `validate_p0400.py` (steps 1–5), `p0400_verify.py` (step 1 alternatives, step 5 cross-check),
`p0400_recording.py` (step 6).

---

## Inputs (Shimadzu datasheet and mount drawing)

| | value |
|---|---|
| grooves | 400 l/mm |
| angle of incidence α | 20.0° |
| entrance arm r | 96.9 mm |
| band | 400–1100 nm, 1st order |
| band edges | r′₁ = 97.1 mm, β₁ = −10.49° (400 nm); r′₂ = 99.5 mm, β₂ = +5.62° (1100 nm) |
| detector | r′ᵢ = 96.9 mm, βᵢ = −7.51° (perpendicular distance and angle of the detector plane) |
| detector length / dispersion | 27.7 mm / 25.3 nm/mm |
| blaze / efficiency | 600 nm nominal, ≥ 50 % relative efficiency there |
| blank | 40 × 40 × 9.5 mm (datasheet); 13.6 mm thick at the corner (mount drawing) |

## Step 1 — the 400 l/mm is at the pole (ν₀ = 400)

Grating equation at the pole, first order: sin β = ν₀ λ − sin α.

| | pole density | β(400 nm) | β(1100 nm) | detector length |
|---|---|---|---|---|
| **400 at the pole** | 400.0 | **−10.49°** | **+5.62°** | **27.66 mm** |
| 400 at the blank edge y = −20 mm | 411.9 | −10.21° | +6.38° | 28.53 mm |
| 400 at the blank edge y = +20 mm | 401.8 | −10.45° | +5.74° | 27.78 mm |
| datasheet | | −10.49° | +5.62° | 27.7 mm |

A change of 1 l/mm at the pole moves β(1100 nm) by 0.063°, so the datasheet angles pin ν₀ = 400 ± 0.16 l/mm.

## Step 2 — the detector is one flat, tilted plane

r′₁ and r′₂ are the two ends of the detector, not a focal curve. The datasheet gives the plane directly:
perpendicular distance r′ᵢ from the grating centre, normal at βᵢ = **−7.51°**, and the distance along any
wavelength's ray is r′(λ) = r′ᵢ / cos(βᵢ − β) (brochure eq. 23).

We use βᵢ = −7.51° as printed and **r′ᵢ = 96.93 mm**. The printed 96.9 is rounded: with βᵢ = −7.51°, the printed
r′₁ = 97.1 and r′₂ = 99.5 are reproduced only for r′ᵢ between 96.92 and 96.95 mm (96.93 gives 97.06 and 99.53).
So the four printed mounting numbers are mutually consistent, and the plane is fully determined by them.

(An earlier version fitted the plane through the rounded r′₁ and r′₂ alone and got 96.96 mm at −7.36°. That
tilt is 0.15° off the printed βᵢ — inside the rounding of r′₁/r′₂ but not of βᵢ. Switching to the printed βᵢ
moved C₁ by 0.4 % and the resolution by < 0.01 nm mean.) The dispersion along this tilted plane is brochure eq. 24,
D = cos β · cos(βᵢ − β) / (N m r′) — 25.28 nm/mm at 400 nm, 24.37 nm/mm at 1100 nm, 25.3 nm/mm on average.

## Step 3 — the radius of curvature from the blank

Not printed, but the blank gives it: sag s = 13.6 − 9.5 = 4.1 mm at the corner, h = 20√2 = 28.28 mm, so

$$R = \frac{h^2 + s^2}{2s} = 99.61\ \text{mm}.$$

(Reading 13.6 mm at the 20 mm edge-midpoint would give 50.8 mm, which cannot focus near the ~97 mm arms.)
Independent cross-checks: freeing R in the ray-trace fit gives 100.00 mm; the classical focus condition gives
100.29 mm — three estimates within 0.7 %.

## Step 4 — what each coefficient does

Along the grating, the number of grooves from the pole is N(y) = ν₀y + C₁y²/2 + C₂y³/3 + C₃y⁴/4. In the
Noda–Namioka light-path function the groove term adds mλN(y) to the optical path, so each coefficient
lands on one aberration:

- **C₁ → defocus.** The defocus coefficient becomes
  $$F_{20}(\lambda) = \frac{\cos^2\alpha}{r} - \frac{\cos\alpha}{R} + \frac{\cos^2\beta}{r'} - \frac{\cos\beta}{R} + \lambda C_1$$
  so defocus vanishes at wavelength λ when **C₁ = −F₂₀,geom(λ) / λ**.
- **C₂ → coma** (adds a term ∝ λC₂ to the coma coefficient F₃₀).
- **C₃ → the next order** (∝ λC₃); it turns out not to be needed here.

(Sign convention: *y* toward the slit side, as in the ray tracer; the analytic formulas' coordinate runs the
other way, which flips the sign of C₁ and C₃ — checked numerically below.)

**The key constraint.** On the fixed datasheet plane with R = 99.61 mm, the C₁ that would cancel defocus is
different at every wavelength:

| λ (nm) | 400 | 500 | 600 | 700 | 800 | 900 | 1000 | 1100 |
|---|---|---|---|---|---|---|---|---|
| C₁ needed (l/mm²) | 0.577 | 0.300 | 0.171 | 0.126 | 0.134 | 0.177 | 0.244 | 0.328 |

A real grating has one C₁. So a single groove pattern — whether a VLS polynomial or a two-point holographic
recording, whose defocus contribution is also λ × constant — can cancel defocus at only about **two
wavelengths**. That is the flat-field compromise, and it is why resolution reaches the floor at two
wavelengths and rises between and beyond them.

## Step 5 — solving C₁, C₂, C₃

With ν₀, the detector plane and R all fixed from the datasheet, the only unknowns are C₁–C₃. They are found
by ray tracing: rays from the slit through the exact spherical surface, diffracted with the local groove
density ν(y), onto the fixed detector plane, at 8 wavelengths (400–1100 nm); C₁–C₃ are adjusted to minimise
the spot size in the dispersion direction (mean + 0.6 × worst).

| | C₁ (l/mm²) | C₂ (l/mm³) | C₃ (l/mm⁴) |
|---|---|---|---|
| fit at f/5.3 | +0.2530 | −0.01715 | ≈ 0 |
| fit at f/10 | +0.2517 | −0.01733 | ≈ 0 |

Defocus nulls fall near **510 nm and 1010 nm**.

**Cross-check (two independent methods).** Fitting C₁ alone to focus one wavelength in the ray tracer
(at f/20, so uncorrected coma does not pull the best focus) reproduces the analytic defocus condition:

| λ (nm) | 400 | 500 | 600 | 700 | 800 | 900 | 1000 | 1100 |
|---|---|---|---|---|---|---|---|---|
| ray tracer | 0.576 | 0.300 | 0.170 | 0.125 | 0.133 | 0.177 | 0.244 | 0.328 |
| analytic | 0.577 | 0.300 | 0.171 | 0.126 | 0.134 | 0.177 | 0.244 | 0.328 |

**The resulting profile** (f/5.3 fit): 388.1 l/mm at y = −20 mm → **400.0 at the pole** → peak **400.9 at
y = +7.4 mm** → 398.2 at y = +20 mm. A chirp of about 3 % across the blank.

**What the fitted grating reproduces:** band-edge angles exact and detector length 27.66 mm (27.7). (The detector
plane is an input taken from the datasheet, so it is not a check.) These test ν₀, not C₁–C₃: where a wavelength
lands is set by the pole density, and the chirp only focuses the rays there. Traced, the centroid detector length
is 27.686 mm with constant grooves and 27.657 mm with the solved chirp at f/5.3 (27.665 / 27.657 at f/10) — both
27.7 mm on the datasheet. So the datasheet cannot distinguish C₁–C₃; Step 7 can. With the relay and a 25 µm slit, the resolution (FWHM of the traced line) is
0.63–0.68 nm at 1000 and 500 nm and 1.43 nm averaged over 400–1100 nm.

## Step 6 — can a real holographic exposure write this?

For seven common recording lasers (266–488 nm), a two-point exposure (two point sources, the pole density
pinned to 400 l/mm) was searched to minimise the band resolution. Every one reproduces the fitted grating:
2.01–2.02 nm mean against 2.02 nm for the free polynomial, groove profiles within 0.1 l/mm across the
illuminated aperture. So the solved chirp is physically writable — but the **recording wavelength cannot be
identified** from optical data: for any λ₀ there is a source layout that writes the same grooves. Solutions
fall into two families (sources on either side of the grating normal) that differ by ~1.3 l/mm at the blank
edges, outside the illuminated aperture.

## Step 7 — the one real unknown, and how to measure it

Everything above is reconstructed from published numbers; the actual grooves are Shimadzu's. On a sample:

1. **Groove map.** Scan a ~1 mm laser beam across the blank at fixed incidence and record the first-order spot
   on a camera 0.5–1 m away; ν(y) = (sin α + sin β)/λ locally. Record the zero-order spot too, to remove the
   local surface tilt of the concave blank. At 633 nm, 1 l/mm moves the spot ~0.3 mm at 500 mm, so C₁ and C₂
   (and which way round the blank goes) are directly measurable.
2. **Line-lamp test.** Hg/Ar lines across the band: the resolution should bottom out near 500 and 1000 nm and
   rise between — the fingerprint of this chirp.
3. **Orientation.** The chirp is asymmetric; mounting the blank rotated 180° takes the model from 2.02 nm to
   ~10 nm mean. Shimadzu's brochure shows a notch on the top side of the blank; confirm which end of the chirp it marks.

## Cross-check against Shimadzu's grating brochure

Shimadzu's brochure page on aberration-corrected flat-field gratings (Fig. 28, eqs. 22–24) gives the mount
geometry but nothing on the groove pattern. What it does confirm:

| Brochure | What it settles here |
|---|---|
| "N: groove density **at the center** of the concave grating" | Step 1: ν₀ = 400 l/mm at the pole — now stated by Shimadzu, not only inferred. |
| eq. 22, β = sin⁻¹(Nmλ − sin α) | Same grating equation used everywhere here. |
| eq. 23, r′ = r′ᵢ / cos(βᵢ − β); r′ᵢ measured perpendicular to the detector | Step 2: one flat, tilted detector plane, now taken directly from the table's r′ᵢ and βᵢ. |
| eq. 24, D = cos β · cos(βᵢ − β) / (N m r′); the quoted dispersion is range / L | Exposed a plate-scale error in our analytic scripts (corrected — see `VALIDATION.md`). Up to 5 % at 1100 nm; ray-traced results unaffected. |
| Fig. 28, "notch (top side)" on the blank | The blank carries an orientation mark. Which end of the chirp sits on the notch side is still not published — Step 7, item 3. |
| "outside the wavelength range, aberrations increase" | Consistent with a chirp optimised for exactly 400–1100 nm. |

Nothing in the brochure changes C₁–C₃: they remain a reconstruction from the mount, to be confirmed by
measurement (Step 7).

