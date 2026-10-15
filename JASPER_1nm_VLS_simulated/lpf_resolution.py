#!/usr/bin/env python3
"""Concave-grating resolution for ALL JASPER designs via the light-path-function (LPF) aberration
theory (Noda/Namioka/Seya 1974) -- the standard analytic method, adopted after Optiland's
sequential ray trace was found unreliable for these off-axis reflective diffractive mounts
(it could match dispersion OR focus but not both, so its spot sizes were fitting artifacts).

LPF is self-consistent: the grating equation sets dispersion; F20 (defocus), F30 (coma) and
F02 (astigmatism) set spot sizes straight from the cited geometry. For each design we report:
  - CLASSICAL grooves (no aberration correction): rigorous, no unknowns -> a resolution floor
  - IDEALISED BOUND: F20 (defocus) nulled at EVERY wavelength, coma-limited. Physically this is a
    detector that follows the curved focal surface; a single groove chirp (VLS or two-point hologram) on a
    FLAT detector nulls defocus at only ~2 wavelengths (see validate_models.py for the real single-chirp case)
and the astigmatic line height (throughput). Only the best-case assumes an ideal groove
correction (unpublished for the catalog gratings); the classical case has no unknowns."""
import numpy as np, json

SLIT=20.0; PIX=8.0  # um

def analyse(name, N_lpm, alpha_deg, rA, rB, Rtan, rho, EPD, band, beta_ref_deg):
    d_nm=1e6/N_lpm; a=np.radians(alpha_deg); bref=np.radians(beta_ref_deg); m=1
    lams=np.array(band)  # um
    def beta(l): return np.arcsin(m*(l*1000)/d_nm - np.sin(a))
    def rp(b):  return rB/np.cos(b-bref)                 # flat detector faces the beta_ref ray
    def plate(b,r): return d_nm*np.cos(b)*np.cos(b-bref)/(m*r)   # nm/mm along the flat detector, Shimadzu eq (24)
    def F20(b,r,R): return np.cos(a)**2/rA-np.cos(a)/R + np.cos(b)**2/r-np.cos(b)/R
    def F30(b,r,R): return (np.sin(a)/rA)*(np.cos(a)**2/rA-np.cos(a)/R)+(np.sin(b)/r)*(np.cos(b)**2/r-np.cos(b)/R)
    def F02(b,r,rho):
        inv=0.0 if not np.isfinite(rho) else (np.cos(a)+np.cos(b))/rho
        return 1.0/rA+1.0/r-inv
    def disp_rms_classical(b,r):
        w=np.linspace(-EPD/2,EPD/2,101)
        dy=(r/np.cos(b))*(F20(b,r,Rtan)*w+1.5*F30(b,r,Rtan)*w**2); dy-=dy.mean()
        return np.sqrt(np.mean(dy**2))*1000
    def disp_rms_corrected(b,r):   # F20 nulled -> coma only
        span=(r/np.cos(b))*1.5*abs(F30(b,r,Rtan))*(EPD/2)**2
        return span*1000/np.sqrt(12)
    out={"name":name,"N":N_lpm,"alpha":alpha_deg,"rA":rA,"rB":rB,"Rtan":Rtan,"rho":rho,"EPD":EPD,"fnum":rA/EPD,"rows":[]}
    for l in lams:
        b=beta(l); r=rp(b); P=plate(b,r); fout=r/EPD; diff=1.028*l*fout
        rc=disp_rms_classical(b,r); rk=disp_rms_corrected(b,r)
        ah=abs(r*F02(b,r,rho)*(EPD/2))*2*1000/1000.0  # mm full height
        dl_c=np.hypot(np.hypot(SLIT,2.355*rc),np.hypot(PIX,diff))/1000*P
        dl_k=np.hypot(np.hypot(SLIT,2.355*rk),np.hypot(PIX,diff))/1000*P
        out["rows"].append(dict(lam=int(l*1000),beta=float(np.degrees(b)),rp=float(r),plate=float(P),
                                rms_classical=float(rc),rms_corr=float(rk),astig_mm=float(ah),
                                dlam_classical=float(dl_c),dlam_corr=float(dl_k)))
    R=out["rows"]
    out["plate_lo"],out["plate_hi"]=R[0]["plate"],R[-1]["plate"]
    out["dlam_corr_mean"]=float(np.mean([x["dlam_corr"] for x in R]))
    out["dlam_corr_worst"]=float(max(x["dlam_corr"] for x in R))
    out["dlam_classical_mean"]=float(np.mean([x["dlam_classical"] for x in R]))
    out["astig_mm"]=float(np.mean([x["astig_mm"] for x in R]))
    return out

def show(o):
    print(f"\n=== {o['name']}  |  {o['N']:.0f} l/mm  alpha={o['alpha']}deg  arms {o['rA']:.1f}/{o['rB']:.1f}  Rtan={o['Rtan']:.1f} rho={o['rho'] if np.isfinite(o['rho']) else 'flat'}  EPD={o['EPD']} (f/{o['fnum']:.1f}) ===")
    print(f"  {'lam':>5}{'beta':>8}{'plate':>7}{'RMScla':>8}{'RMScorr':>8}{'astig(mm)':>10}{'dlamCorr':>9}")
    for x in o["rows"]:
        print(f"  {x['lam']:>5}{x['beta']:>8.2f}{x['plate']:>7.1f}{x['rms_classical']:>8.1f}{x['rms_corr']:>8.1f}{x['astig_mm']:>10.2f}{x['dlam_corr']:>9.2f}")
    print(f"  >> plate {o['plate_lo']:.1f}-{o['plate_hi']:.1f} nm/mm | idealised-bound MEAN {o['dlam_corr_mean']:.2f} WORST {o['dlam_corr_worst']:.2f} nm | classical MEAN {o['dlam_classical_mean']:.2f} nm | astig line {o['astig_mm']:.2f} mm")

if __name__=="__main__":
    band=[0.5,0.6,0.7,0.8,0.9,1.0]
    # --- Simulated concave VLS grating (our theoretical design), sag radius for the toroid nulls astig ---
    a=np.radians(12.34); rAs=137.357; rBs=126.75; Rs=133.0
    bc=np.arcsin((0.75*1000)/(1e6/230.0)-np.sin(a)); rho_sim=(np.cos(a)+np.cos(bc))/(1.0/rAs+1.0/rBs)
    designs=[
      analyse("Simulated VLS — SPHERICAL",  230, 12.34, 137.357,126.75, 133.0, 133.0,   26.0, band, -2.4),
      analyse("Simulated VLS — TOROIDAL",   230, 12.34, 137.357,126.75, 133.0, rho_sim, 26.0, band, -2.4),
      analyse("Shimadzu P0400-01 (real, spherical)", 400, 20.0, 96.9, 96.93, 99.61, 99.61, 18.3, band, -7.51),  # R from the blank sag, plane from the datasheet r'i / beta_i (validate_p0400.py)
    ]
    print(f"(simulated toroid sag radius rho={rho_sim:.1f} mm)")
    for o in designs: show(o)
    json.dump(designs, open("lpf_all_designs.json","w"), indent=1)
    print("\nsaved lpf_all_designs.json")
