"""Glass dispersion + catalog achromat prescriptions + a real-surface 3-D ray tracer for the fiber->slit relay."""
import numpy as np
# Schott Sellmeier coefficients (B1,B2,B3,C1,C2,C3), lambda in um
GLASS = {
 "N-BK7":  (1.03961212, 0.231792344, 1.01046945, 0.00600069867, 0.0200179144, 103.560653),
 "N-SF5":  (1.52481889, 0.187085527, 1.42729015, 0.011254756,   0.0588995392, 129.141675),
 "N-BAF10":(1.5851495,  0.143559385, 1.08521269, 0.00926681282, 0.0424489805, 105.613573),
 "N-SF10": (1.62153902, 0.256287842, 1.64447552, 0.0122241457,  0.0595736775, 147.468793),
}
CATALOG_ND_VD = {"N-BK7": (1.51680, 64.17), "N-SF5": (1.67270, 32.25), "N-BAF10": (1.67003, 47.11), "N-SF10": (1.72828, 28.53)}
def n(glass, lam_um):
    if glass == "air": return 1.0
    B1, B2, B3, C1, C2, C3 = GLASS[glass]; l2 = lam_um ** 2
    return np.sqrt(1 + B1*l2/(l2-C1) + B2*l2/(l2-C2) + B3*l2/(l2-C3))

# Edmund TECHSPEC VIS-NIR achromats (R1 faces the long / infinite conjugate). Surfaces: (radius, thickness after, glass after)
LENS = {
 # 25 mm diameter (clear aperture 24 mm)
 "49-352": dict(efl=30.0,  bfl=22.23,  R=(21.17, -16.08, -118.66), t=(11.04, 3.00), g=("N-BAF10", "N-SF10"), dia=25.0, ca=24.0, inr=13923),
 "49-361": dict(efl=125.0, bfl=120.89, R=(76.28, -54.55, -162.43), t=(6.00, 2.40),  g=("N-BK7", "N-SF5"),   dia=25.0, ca=24.0, inr=13721),
 "49-360": dict(efl=100.0, bfl=95.92,  R=(61.47, -44.64, -129.94), t=(6.00, 2.50),  g=("N-BK7", "N-SF5"),   dia=25.0, ca=24.0, inr=None),
 # 12.5 mm diameter (clear aperture 11.5 mm)
 "49-321": dict(efl=14.0,  bfl=9.92,   R=(10.60, -7.72, -56.54),   t=(6.23, 1.50),  g=("N-BAF10", "N-SF5"),  dia=12.5, ca=11.5, inr=9736),
 "49-329": dict(efl=60.0,  bfl=57.59,  R=(36.51, -26.33, -78.35),  t=(3.50, 1.50),  g=("N-BK7", "N-SF5"),    dia=12.5, ca=11.5, inr=9736),
}
def surfaces(lens, reverse=False):
    """list of (radius, thickness_to_next, glass_after) in the order light meets them."""
    L = LENS[lens]; R1, R2, R3 = L["R"]; t1, t2 = L["t"]; g1, g2 = L["g"]
    if not reverse: return [(R1, t1, g1), (R2, t2, g2), (R3, 0.0, "air")]
    return [(-R3, t2, g2), (-R2, t1, g1), (-R1, 0.0, "air")]

def paraxial(surfs, lam_um):
    """EFL and BFL of a surface list (object at infinity) by y-nu paraxial trace."""
    y, u, nprev = 1.0, 0.0, 1.0
    for i, (Rr, t, g) in enumerate(surfs):
        n2 = n(g, lam_um); phi = (n2 - nprev) / Rr
        u = (nprev * u - y * phi) / n2
        if i < len(surfs) - 1: y = y + t * u
        nprev = n2
    return -1.0 / u, -y / u      # EFL, BFL (from last vertex)

def trace(P, D, surfs, z0, lam_um, ca_r=12.0):
    """Real-ray trace. P,D: (N,3) start points / unit directions (z = optical axis). Surfaces start at z0."""
    P = P.copy(); D = D.copy(); z = z0; nprev = 1.0; alive = np.ones(len(P), bool)
    for Rr, t, g in surfs:
        C = np.array([0, 0, z + Rr])                         # sphere centre
        oc = P - C; b = np.sum(D * oc, 1); c = np.sum(oc * oc, 1) - Rr * Rr; disc = b * b - c
        alive &= disc >= 0
        s = -b - np.sign(Rr) * np.sqrt(np.clip(disc, 0, None))   # nearer intersection for the vertex cap
        P = P + s[:, None] * D
        N = (P - C) / Rr; N = -N if Rr > 0 else N               # normal pointing back toward the incoming light... make it face -D
        N = np.where(np.sum(N * D, 1)[:, None] > 0, -N, N)
        n2 = n(g, lam_um); mu = nprev / n2; ci = -np.sum(N * D, 1)
        k = 1 - mu * mu * (1 - ci * ci); alive &= k >= 0
        D = mu * D + (mu * ci - np.sqrt(np.clip(k, 0, None)))[:, None] * N
        D /= np.linalg.norm(D, axis=1, keepdims=True)
        alive &= np.hypot(P[:, 0], P[:, 1]) <= ca_r            # clear-aperture radius
        z += t; nprev = n2
    return P, D, alive, z

if __name__ == "__main__":
    print("glass check (Sellmeier vs catalog):")
    for g in GLASS:
        nd = n(g, 0.5875618); nF = n(g, 0.4861327); nC = n(g, 0.6562725); vd = (nd - 1) / (nF - nC)
        print(f"  {g:8s} nd {nd:.5f} (cat {CATALOG_ND_VD[g][0]:.5f})  vd {vd:5.2f} (cat {CATALOG_ND_VD[g][1]:.2f})")
    print("lens check (paraxial @ 587.6 nm vs Edmund datasheet):")
    for L in LENS:
        efl, bfl = paraxial(surfaces(L), 0.5876)
        print(f"  #{L}: EFL {efl:7.2f} (spec {LENS[L]['efl']:.2f})  BFL {bfl:7.2f} (spec {LENS[L]['bfl']:.2f})")
