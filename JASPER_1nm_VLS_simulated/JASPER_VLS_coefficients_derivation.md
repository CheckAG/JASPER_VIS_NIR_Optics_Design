# Deriving the VLS groove coefficients from a real holographic recording

*JASPER optics — supplement to the simulated ~1 nm concave VLS grating.*
*How the polynomial groove chirp ν(y) = 230 + C₁y + C₂y² + C₃y³ comes from a physically buildable
two-beam exposure, with the equations and a worked example.*

---

## 1. Why this note exists

The main model corrects field aberration with a variable groove density written as a polynomial

$$\nu(y) = \nu_0 + C_1 y + C_2 y^2 + C_3 y^3 \quad [\text{lines/mm}]$$

and fits C₁, C₂, C₃ to minimise the spot size. That is convenient, but a real grating is **not** made
by dialling three independent numbers. It is made by exposing the blank to **two coherent point
sources** C and D at a recording wavelength λ₀. Those five numbers — the source distances and angles
(r_C, γ, r_D, δ) and λ₀ — fix the whole groove function. So C₁, C₂, C₃ are **not free**; they are
locked to each other by the recording. This note derives that link, so the "best case" is a target a
real exposure can actually hit, not an unconstrained curve.

---

## 2. The exact groove function

The grooves of a holographic grating are the interference fringes of the two sources. The local groove
density at a surface point **P** is

$$\boxed{\;\nu(\mathbf P) = \frac{1}{\lambda_0}\,\big(\hat{\mathbf u}_C - \hat{\mathbf u}_D\big)\cdot \hat{\mathbf g}(\mathbf P)\;}$$

where

- $\hat{\mathbf u}_C = (\mathbf P-\mathbf C)/|\mathbf P-\mathbf C|$ and $\hat{\mathbf u}_D = (\mathbf P-\mathbf D)/|\mathbf P-\mathbf D|$ are the unit vectors from each source to the point, and
- $\hat{\mathbf g}(\mathbf P)$ is the in‑surface unit vector perpendicular to the grooves (the dispersion tangent).

This is the Noda–Namioka result; it is exact for any source positions and any substrate shape. Every
number below follows from it.

Geometry (dispersion plane y–z, grooves along x, pole at the origin, normal +z): a source at distance
r and angle θ from the normal is at $\mathbf C=(0,\;r\sin\theta,\;r\cos\theta)$.

---

## 3. The mean line density — set by the angles alone

At the pole **P** = 0 the tangent is $\hat{\mathbf g}=\hat{\mathbf y}$ and $\hat{\mathbf u}_C=-\mathbf C/r_C=(0,-\sin\gamma,-\cos\gamma)$, so

$$\boxed{\;\nu_0 = \frac{\sin\delta - \sin\gamma}{\lambda_0}\;}$$

The **mean groove density is fixed purely by the two source angles**, not their distances. For the
JASPER design this is the constraint that pins 230 l/mm:

$$\frac{\sin 31.25^\circ - \sin 23.98^\circ}{488\ \text{nm}} = \frac{0.5189-0.4066}{488\times10^{-6}\,\text{mm}} = 230.0\ \text{l/mm}.$$

This equation removes one degree of freedom: once you pick λ₀ and one angle, the other angle is set.

---

## 4. The chirp — the plane-grating skeleton

Expand the source-to-point distance for a point a distance *w* along a **flat** grating
($\mathbf P=(0,w,0)$). With $|\mathbf P-\mathbf C| = \sqrt{w^2 - 2w\,r_C\sin\gamma + r_C^2}$ and the
series $\sqrt{1+b\varepsilon+\varepsilon^2}=1+\tfrac{b}{2}\varepsilon+\tfrac{4-b^2}{8}\varepsilon^2+\dots$,

$$|\mathbf P-\mathbf C| = r_C - \sin\gamma\,w + \frac{\cos^2\gamma}{2 r_C}\,w^2 + \frac{\sin\gamma\cos^2\gamma}{2 r_C^2}\,w^3 + \dots$$

The groove count from the pole is $N(w) = \tfrac{1}{\lambda_0}\big[(|\mathbf P-\mathbf C|-r_C) - (|\mathbf P-\mathbf D|-r_D)\big]$,
and the density is $\nu(w)=\mathrm dN/\mathrm dw$. Collecting terms:

$$\boxed{\;C_1^{\text{plane}} = \frac{1}{\lambda_0}\!\left(\frac{\cos^2\gamma}{r_C}-\frac{\cos^2\delta}{r_D}\right),\qquad
C_2^{\text{plane}} = \frac{3}{2\lambda_0}\!\left(\frac{\sin\gamma\cos^2\gamma}{r_C^2}-\frac{\sin\delta\cos^2\delta}{r_D^2}\right)\;}$$

So the **chirp is driven by the source distances** (r_C, r_D), while the mean density was driven by the
angles. C₃ follows from the next term in the same series. These formulas are exact for a plane
holographic grating — verified to the digit against a direct numerical evaluation of §2.

---

## 5. The concave correction — and why it can't be skipped

On a **concave** blank the point P has sag, the surface normal rotates as you move along it, and the
groove tangent $\hat{\mathbf g}$ tilts with it. Each of those adds terms linear and quadratic in *w*.
For the JASPER substrate (R = 133 mm, aperture ±13 mm) these curvature terms **dominate** the chirp —
they even flip C₁'s sign:

| coefficient | plane-grating formula (§4) | concave R = 133 mm (delivered) |
|---|---:|---:|
| ν₀ (l/mm)     | 230.000  | 230.000 |
| C₁ (l·mm⁻²)   | +0.813   | **−0.092** |
| C₂ (l·mm⁻³)   | −3.6e‑4  | −2.9e‑3 |
| C₃ (l·mm⁻⁴)   | −1.1e‑5  | −6e‑6 |

(Both columns are for the same recording: λ₀ = 488 nm, C = (r 329.2 mm, +23.98°), D = (r 341.7 mm,
+31.25°).) The lesson: **ν₀ is exactly the angle formula on any substrate, but the chirp on a concave
grating must be taken from the exact expression of §2 evaluated on the real surface** — a five-line
computation, not the plane formula. That is exactly what `recording_optimize_sim.py` / the `nu_holo`
function do; the polynomial C₁, C₂, C₃ are then just the Taylor coefficients of that exact ν(y) about
the pole.

---

## 6. Forward and inverse

- **Forward** (what a given exposure produces): pick (r_C, γ, r_D, δ, λ₀) → evaluate §2 on the surface →
  read off ν₀, C₁, C₂, C₃.
- **Inverse** (how we design): the pole constraint $\sin\delta = \sin\gamma + \nu_0\lambda_0$ removes one
  DOF; search the remaining four over **real recording-laser lines** to minimise the band resolution
  directly. A free polynomial can name any C₁, C₂, C₃; a two‑beam exposure can only reach the subset
  spanned by (r_C, γ, r_D, δ, λ₀). The optimisation confirms the design's ~0.9 nm target lies inside
  that reachable subset — so it is buildable, not just a curve fit.

## 7. Limits (stated honestly)

A simple two‑point spherical‑wave exposure buys only **low‑order** correction (essentially the C₁, C₂
seen here). If a design needs more — strong higher‑order coma across a very wide field — the recording
itself must be made with an **aspheric wavefront, an auxiliary optic, or a computer‑generated
hologram**. Those are more capable and more expensive, and are a deliberate fork in the design, not a
free upgrade. For JASPER's 500–1000 nm band at f/10 the simple two‑beam recording is enough.

---

*Reproduce: `python recording_optimize_sim.py` prints the recording geometry, the ν(y) it produces, and
the resolution it reaches, and writes `jasper_recording_optimize.png` (physical vs free‑polynomial ν(y)
and the per‑band resolution) and `jasper_recording_spots.png` (the ray‑traced spots for the real
recording). The coefficient table above is `recording_optimize_sim.json` Taylor‑expanded about the pole.*
