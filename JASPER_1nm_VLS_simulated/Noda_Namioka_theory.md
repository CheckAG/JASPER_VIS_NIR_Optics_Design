# The Noda–Namioka Light-Path-Function Theory of Concave Gratings

*A working write-up for the JASPER spectrometer design — the analytic method behind all the
resolution and astigmatism numbers in the design record.*

CheckAg Private Limited · Tony Francis

---

## 1. Why this theory

A concave diffraction grating does two jobs at once: it **disperses** light (splits it by
wavelength) and **focuses** it (forms an image of the slit on the detector), with no separate
mirrors. Because dispersion and focusing happen at the *same* curved, ruled surface, you cannot
treat them separately — you need a theory that handles the grating equation and the imaging
aberrations together, self-consistently.

That theory is the **light-path function** (LPF) formalism of Noda, Namioka and Seya (1974). It is
the standard tool every grating manufacturer (Shimadzu, HORIBA, Newport) uses to design
aberration-corrected concave gratings, and it is what we used for JASPER after the sequential ray
trace proved unreliable for these off-axis reflective mounts.

Its great virtue: dispersion falls out of one term (the grating equation), and every aberration —
defocus, coma, astigmatism — falls out of the *other* terms of a single power series, **computed
directly from the mounting geometry**. No ray trace, no free parameters.

---

## 2. The light-path function

Set up coordinates with the grating pole (centre) at the origin **O**, the grating normal along
*x*, the dispersion plane = the *xy*-plane, and the grooves running parallel to *z* at the pole.

- Entrance slit point **A**: distance *r* from O, at angle *α* from the normal.
- Image point **B**: distance *r′*, at angle *β* from the normal.
- A general point **P** on the grating surface has in-plane coordinates **w** (along *y*, the
  dispersion direction) and **l** (along *z*, the groove/sagittal direction). The surface sag bends
  P off the tangent plane by an amount set by the substrate radii.

By Fermat's principle, the actual ray from A to B via P makes the **optical path length
stationary**. For a grating, the path also picks up a phase term from the grooves. The
**light-path function** is:

```
F(w, l)  =  <AP>  +  <PB>  +  (m·λ / d) · n(w)
```

where `<AP>`, `<PB>` are the geometric distances, *m* is the diffraction order, *λ* the wavelength,
*d* the groove spacing at the pole, and `n(w)` counts the grooves crossed out to position *w*. (For
a **holographic** grating the groove term is written instead as `(m·λ/λ₀)·(<CP> − <DP>)`, where
C and D are the two recording point sources and λ₀ the recording wavelength — this is what lets the
grooves be curved and unequally spaced.)

The imaged ray positions are found by expanding F in powers of the aperture coordinates *w* and *l*
and requiring it to be stationary.

---

## 3. The power series and the aberration coefficients

Expand F as a power series in the grating-aperture coordinates:

```
F(w,l) = F₀₀ + F₁₀·w + ½F₂₀·w² + ½F₀₂·l² + ½F₃₀·w³ + ½F₁₂·w·l² + …
```

Each coefficient F_{ij} has a direct physical meaning. For a **toroidal** grating with meridional
(tangential) radius *R* and sagittal radius *ρ*, classically-ruled grooves of density 1/*d* at the
pole:

**F₁₀ — the grating equation (dispersion).**
```
F₁₀ = sin α + sin β − m·λ/d
```
Setting **F₁₀ = 0** gives the grating equation `sin α + sin β = m·λ/d`, which fixes the diffraction
angle β for each wavelength. This term *is* the dispersion; everything below is imaging.

**F₂₀ — defocus (meridional / spectral focusing).**
```
F₂₀ = (cos²α/r − cos α/R) + (cos²β/r′ − cos β/R)
```
Setting **F₂₀ = 0** is the **tangential focus condition** — it places the sharp spectral focus on
the detector. A plain spherical grating can satisfy it at only *one* wavelength; the whole art of a
**flat-field** grating is keeping F₂₀ ≈ 0 across the *whole* band (see §5).

**F₀₂ — astigmatism (sagittal focusing).**
```
F₀₂ = (1/r − cos α/ρ) + (1/r′ − cos β/ρ)
```
Setting **F₀₂ = 0** is the **sagittal focus condition**. A spherical grating (ρ = R) generally
*cannot* satisfy F₂₀ = 0 and F₀₂ = 0 at the same time — the leftover is **astigmatism**, which
stretches each wavelength into a vertical line. A **toroidal** grating uses a separate ρ to null
F₀₂ and remove it. (This is the entire spherical-vs-toroidal throughput story in the design record.)

*Solving for the ideal toroid ρ.* Set F₀₂ = 0 and rearrange for ρ:
```
(cos α + cos β)/ρ = 1/r + 1/r′        →        ρ = (cos α + cos β) / (1/r + 1/r′)
```
This is a closed-form solution — no optimizer needed. But β depends on wavelength, so F₀₂ can only
be nulled *exactly* at one wavelength; ρ is chosen at the band centre, leaving a small residual
astigmatism at the band edges. That is exactly what `lpf_resolution.py` computes: it evaluates β at
the mid-band wavelength and sets `rho = (cosα + cosβ_centre)/(1/rA + 1/rB)`. (For JASPER's simulated
230 l/mm design this gives ρ ≈ 130 mm against a tangential R = 133 mm — i.e. the toroid is only ~2%
"out of round," which is precisely why it's optically effective yet mechanically hard to make, §5.)

**F₁₂ — astigmatic coma (the w·l² term).**
```
F₁₂ = (sin α/r)(1/r − cos α/ρ) + (sin β/r′)(1/r′ − cos β/ρ)
```
This couples the two planes: it makes the astigmatic line *curved* and makes the astigmatism vary
across the dispersion aperture. It is the leading residual once F₂₀ and F₀₂ are corrected, alongside
F₃₀. For JASPER's near-normal mounts it is small, but a full holographic design minimizes it too.

**F₃₀ — primary coma.**
```
F₃₀ = (sin α/r)(cos²α/r − cos α/R) + (sin β/r′)(cos²β/r′ − cos β/R)
```
Coma is the dominant *residual* aberration once defocus is corrected. Crucially, **F₃₀ depends only
on the mounting geometry (r, r′, α, β, R)** — not on the groove-correction profile — which is why
the coma-limited resolution is defensible even when a vendor doesn't publish the recording data.

Higher terms (F₄₀ spherical aberration, F₂₂, …) exist but are usually small for these mounts.

---

## 4. From coefficients to spot sizes and resolution

The **transverse ray aberrations** — how far a ray through aperture point (w, l) lands from the
ideal image point, in the detector plane — are the derivatives of F:

```
Δy  (dispersion axis)   =  (r′ / cos β) · ∂F/∂w  =  (r′/cos β)·( F₂₀·w + 3⁄2·F₃₀·w² + … )
Δz  (sagittal axis)     =   r′ · ∂F/∂l           =   r′·( F₀₂·l + … )
```

(The F₁₀ term vanishes at the design β.) Sweeping (w, l) over the illuminated grating aperture
traces out the **spot diagram**:

- The spread of **Δy** is the blur along the dispersion direction → it sets **spectral resolution**.
- The span of **Δz** is the **astigmatic line height** → it sets **throughput** (how much of each
  wavelength's line lands on the detector pixel), *not* resolution.

Converting the dispersion blur to a spectral figure uses the **reciprocal linear dispersion**:
```
plate scale  =  d·cos β / (m·r′)                     [detector perpendicular to the λ ray]
plate scale  =  d·cos β · cos(βᵢ − β) / (m·r′)       [flat detector tilted to the βᵢ ray, nm/mm along it]
```
The second form is Shimadzu's eq. 24 for flat-field gratings, with r′ = r′ᵢ / cos(βᵢ − β) the distance along
the λ ray and r′ᵢ the perpendicular distance to the detector. All flat-detector scripts use it.
and JASPER's resolution budget combines the independent contributions in quadrature:
```
Δλ  =  √( slit² + (2.355·RMS_disp)² + pixel² + FWHM_diffraction² )  ×  plate_scale
```
where RMS_disp is the RMS of Δy over the aperture, and the diffraction term (added separately, since
the LPF is purely geometric) is the Airy/FWHM at the working f-number.

---

## 5. Holographic / VLS gratings — the flat-field trick

A classically-ruled spherical grating satisfies F₂₀ = 0 at only one wavelength, so its resolution
degrades toward the band edges. A **holographic** grating is recorded by interfering two point
sources (C, D), which produces **curved, unequally-spaced grooves** — a variable line spacing (VLS).
Equivalently, the local groove density varies across the aperture as
`ν(w) = ν₀ + ν₁·w + ν₂·w² + …`.

That varying density adds its own contributions to every coefficient — δF₂₀, δF₃₀, δF₀₂ … — with the
δ's set by the recording geometry (r_C, r_D, and the recording angles). The designer chooses the
recording so that:

- **F₂₀ ≈ 0 across the whole band** (flat spectral focus — a "flat-field" grating), and
- F₃₀, F₁₂ are minimized as far as the remaining freedom allows.

This is exactly why the Shimadzu / Lensea polychromator gratings are described as
*"unequally-spaced curved-groove"* holographic gratings. The correction lives in the groove pattern,
which vendors treat as proprietary and rarely publish — hence, when we model a catalog grating, we
can compute the geometry-fixed terms exactly but must treat the correction as *ideal* (best-case).

> **What the code actually does (important).** `lpf_resolution.py` does **not** derive the δF
> recording terms — that algebra (solving for the two recording points C, D that null F₂₀ across the
> band) is the hard part of a full holographic *design*, and we did not implement it. Instead the
> code takes the **best-case shortcut**: it computes the geometry-fixed coefficients F₂₀ (classical,
> uncorrected), F₃₀ (coma) and F₀₂ (astigmatism), and for the "best-case flat-field" number it simply
> *assumes* the correction nulls F₂₀ and reports the coma-limited residual. So there is a deliberate
> gap between this write-up (the full theory) and the code (a two-line best-case/worst-case bracket).
> Running the code, you step through: grating equation → F₂₀/F₃₀/F₀₂ → quadrature — no δF's to wade
> through. Deriving the actual recording geometry (the δF's) is a separate, larger task.

> **Correction (Release 2): the "best-case flat-field" is an idealised bound, not a buildable grating.**
> Nulling F₂₀ at *every* wavelength would need the groove correction itself to change with wavelength.
> It can't: a grating has one fixed groove pattern for all wavelengths, and its defocus contribution is
> λ × (a constant set by the chirp c₁, or by the recording geometry). The geometric defocus on a flat
> detector does not scale like λ, so one chirp can cancel it at only ~2 wavelengths (a third with the
> right detector tilt). The bound is still physically meaningful: it is what a detector *bent to the true
> focal curve* would see, or a scanning instrument that refocuses per wavelength. `validate_models.py`
> reports both — the bound and the real single chirp — in each model, and the recording geometry is now
> derived in `recording_optimize_sim.py` / `JASPER_VLS_coefficients_derivation.md`.

---

## 6. How JASPER uses it

`lpf_resolution.py` implements the coefficients of §3 directly. For each design and wavelength it:

1. solves F₁₀ = 0 for β (dispersion, plate scale);
2. evaluates **F₂₀** for the raw, uncorrected defocus → the *classical* resolution floor;
3. for the **best-case flat-field** number, sets F₂₀ → 0 (ideal groove correction) and takes the
   residual blur from **F₃₀** (coma) — a geometry-only floor, independent of the unpublished groove
   profile;
4. evaluates **F₀₂** for the astigmatic line height (throughput), spherical (ρ = R) vs toroidal
   (ρ nulls F₀₂);
5. folds slit, pixel and diffraction into the quadrature budget above.

Because the grating equation and the aberrations come from *one* self-consistent function, the model
reproduces the Shimadzu P0400-01 published dispersion (25.3 nm/mm) exactly — the check the ray trace
failed.

**Worked example — Shimadzu P0400-01** (400 l/mm, α = 20°, r = 96.9 mm, r′ ≈ 97 mm, R ≈ 100 mm):
F₁₀ = 0 → β = −10.5° (400 nm) to +5.6° (1100 nm), plate 25.3 nm/mm ✓; F₂₀ across 500–1000 nm is
small but non-zero → classical ~1.5 nm; with F₂₀ nulled the F₃₀ coma gives RMS ~3–4 µm → **~0.60 nm**
best-case; F₀₂ on the spherical blank gives a ~2.3 mm astigmatic line.

---

## 7. Assumptions and limitations

- **Geometric only.** The LPF is ray optics; diffraction (the Airy limit) is added separately in
  quadrature. Fine for these ~f/5 mounts where the slit dominates.
- **Power-series truncation.** Valid for moderate apertures; very fast (low f-number) or very wide
  gratings would need higher-order terms.
- **The correction is an input.** The best-case resolution assumes the recording/VLS profile nulls
  F₂₀ ideally. For a catalog grating whose groove data is unpublished, the true resolution lies
  between the classical floor and the best-case; only the geometry-fixed terms (dispersion, coma,
  astigmatism-vs-substrate) are certain.

---

## 8. Independent validation — the standalone ray tracer

The LPF is a *3rd-order* theory. To check it (and to get real ray-geometry spot diagrams), we
built a **standalone geometric ray tracer** — `vls_raytrace.py`, pure NumPy, no Optiland. It
traces thousands of rays from the slit, intersects the **exact** spherical/toroidal grating
surface, applies the **local** variable-line-spacing vector grating equation at each hit point,
and propagates to a flat detector placed on the diffracted chief ray at best focus.

**Why standalone, not Optiland?** Optiland's sequential trace of this off-axis reflective
concave-grating mount was unreliable — for the Shimadzu P0400-01 it could match the published
dispersion *or* the focus but not both at once, so its spot sizes were artifacts. The purpose-built
tracer sidesteps all of Optiland's tilt/coordinate conventions and is transparent end to end.

**What it confirms — two independent methods agree on the hard geometry:**
- Dispersion: **25.3 nm/mm** (ray trace = LPF = Shimadzu datasheet, exactly).
- Astigmatism: the astigmatic line height matches the LPF F₀₂ prediction (~2.3 mm spherical P0400).

**What it corrects — the resolution:** the full ray trace captures *higher-order* aberration
(5th order and up) that the 3rd-order LPF truncates. It is real (an exact sphere and a paraboloid
give the same answer, so it is not a surface artifact) and it is **strongly f-number-dependent**:

| f-number | ray-traced Δλ (mean, 500–1000 nm) |
|---|---|
| f/5.3 | ~2.3 nm |
| f/8   | ~1.35 nm |
| f/10  | ~1.1 nm |
| f/12  | ~1.0 nm |
| (aberration-free floor) | ~0.6 nm |

So the LPF's coma-limited idealised bound (~0.6 nm) is really the **aberration-free floor** — reached
at a single chirp's defocus nulls, and across the band only when the aperture is stopped down. At a usable, bright f/5.3 the grating is
higher-order-aberration-limited (~2 nm). The honest headline is the tradeoff, not a single number:
**sub-nm capable when stopped to ~f/10–12; ~1 nm at f/8; ~2 nm at the bright f/5.3.** f-number is a
real design knob the 3rd-order theory hid. (Note: the ~f/10 setting uses a *smaller* beam on the
grating, not a larger grating — it is set by the aperture stop, not the grating size.)

The tracer therefore does two jobs: it **validates** the analytic model where the analytic model is
exact (dispersion, astigmatism), and it **supersedes** it where the 3rd-order truncation is
optimistic (resolution at fast apertures). Both are kept in the repo; the ray tracer is the
authority for spot size and resolution, the LPF for the closed-form focus/astigmatism conditions.

## 9. References — the original papers

The core theory is the 1974 JOSA trilogy:

1. **H. Noda, T. Namioka, and M. Seya, "Geometric theory of the grating," *J. Opt. Soc. Am.* **64**(8), 1031–1036 (1974).** — the light-path function and the aberration coefficients F_{ij}. *This is the primary reference.*
   https://opg.optica.org/josa/abstract.cfm?uri=josa-64-8-1031
2. H. Noda, T. Namioka, and M. Seya, "Ray tracing through holographic gratings," *J. Opt. Soc. Am.* **64**(8), 1037–1042 (1974). — exact ray-tracing companion.
   https://opg.optica.org/josa/abstract.cfm?uri=josa-64-8-1037
3. H. Noda, T. Namioka, and M. Seya, "Design of holographic concave gratings for Seya–Namioka monochromators," *J. Opt. Soc. Am.* **64**(8), 1043–1048 (1974). — the holographic-recording (VLS) design.
   https://opg.optica.org/josa/abstract.cfm?uri=josa-64-8-1043

Foundational earlier work:

4. T. Namioka, "Theory of the concave grating. I," *J. Opt. Soc. Am.* **49**, 446–460 (1959); and "…III. Seya–Namioka monochromator," *J. Opt. Soc. Am.* **49**, 951–961 (1959).

Modern reference text (uses the same F_{ij} formalism, convenient tables):

5. C. Palmer and E. Loewen, *Diffraction Grating Handbook*, Newport / Richardson Gratings — the "Imaging Properties of Concave Gratings" chapter.

---

*Note on sign conventions:* different texts absorb factors of ½ and signs into the coefficients
differently (Namioka 1959 vs Noda 1974 vs the Grating Handbook). The forms above match the
implementation in `lpf_resolution.py`; when cross-checking against a paper, verify its convention
for the ½-factors and the sign of β.
