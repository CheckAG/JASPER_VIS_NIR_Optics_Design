#!/usr/bin/env python3
"""VALIDATION against a REAL commercial part: the Shimadzu P0400-01 flat-field concave grating.

P0400-01 is a SPHERICAL-substrate, aberration-corrected (holographic / variable-groove) flat-field
grating. (Shimadzu's toroidal flat-field gratings are a separate product line.)

Published (Shimadzu spec table / mount drawing):
  400 grooves/mm, 400-1100 nm, alpha = 20.0 deg, entrance arm r = 96.9 mm,
  exit distances r'1 = 97.1 mm (400 nm, beta1 = -10.49) and r'2 = 99.5 mm (1100 nm, beta2 = +5.62),
  detector r'i = 96.9 mm at beta_i = -7.51 deg, detector length 27.7 mm, dispersion 25.3 nm/mm.
  Radius of curvature R is NOT published.

Detector: Shimadzu's mount (brochure eq. 23) is ONE flat plane, perpendicular distance r'i from the grating
centre along beta_i; the exit distance at any wavelength is r' = r'i / cos(beta_i - beta). We take
beta_i = -7.51 deg straight from the table and r'i = 96.93 mm: the printed 96.9 refined to the only range
(96.92-96.95 mm) that also reproduces the printed r'1 = 97.1 and r'2 = 99.5 after rounding. The detector is
placed there and never moved.

Radius of curvature: not printed, but the blank drawing gives it. The 40 x 40 mm blank is 9.5 mm thick
at the centre and 13.6 mm at the corner (20*sqrt(2) = 28.28 mm out): sag s = 4.1 mm, so
R = (h^2 + s^2) / (2 s) = 99.61 mm. (Reading 13.6 mm at the 20 mm edge-midpoint instead gives 50.8 mm,
which cannot focus at the ~97 mm published arms.) A free ray-trace fit of R lands at ~99.9 mm, and the
classical tangential-focus condition gives 100.29 mm -- three independent estimates within 0.7%.

With R fixed from the blank, the ONLY thing solved is the groove chirp
  nu(y) = 400 + c1*y + c2*y^2 + c3*y^3      (400 l/mm AT THE POLE, y = 0)
by minimising the traced spot size on the fixed datasheet plane across 400-1100 nm.

Physics note: on a fixed flat plane, a single groove chirp (VLS polynomial OR a two-point holographic
recording -- both add a defocus correction proportional to lambda x constant) can null defocus at only
~2 wavelengths. Between/beyond those nulls a residual defocus remains. That is the intrinsic flat-field
compromise, and it is why the resolution is ~floor at two wavelengths and rises elsewhere.
"""
import numpy as np, json
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import minimize, minimize_scalar
import vls_raytrace as R

N=400.0; ALPHA=20.0; rA=96.9
BAND=[0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]; lam_nm=np.array(BAND)*1e3
PUB=dict(L=27.7, r1=(0.4,97.1), r2=(1.1,99.5), disp=25.3, ri=96.9, bi=-7.51)
RI=96.93     # printed r'i = 96.9; 96.92-96.95 also reproduces the printed r'1/r'2 (see docstring)
BLANK=dict(t_centre=9.5, t_corner=13.6, half_side=20.0)          # mm, from the mount drawing
def radius_from_blank(b=BLANK):
    s=b["t_corner"]-b["t_centre"]; h=b["half_side"]*np.sqrt(2)    # corner is 20*sqrt(2) from the centre
    return (h*h+s*s)/(2*s)
RSUB=radius_from_blank()                                           # 99.61 mm
SLIT=20.0; PIX=8.0
a=np.radians(ALPHA)

def beta(l): return np.arcsin(l*1e-3*N-np.sin(a))

def datasheet_plane(ri=RI, bi=PUB["bi"]):
    """Detector plane straight from the datasheet mounting parameters r'i, beta_i (brochure eq. 23)."""
    b=np.radians(bi); n=np.array([0,np.sin(b),np.cos(b)])
    ndet=-n; ex=np.array([1.,0,0]); ed=np.cross(ndet,ex); ed/=np.linalg.norm(ed)
    return (ri*n,ndet,ed,ex), ri, bi

CS, PDIST, PBETA = datasheet_plane()
def plane_r(l): return PDIST/np.cos(beta(l)-np.radians(PBETA))    # distance pole->plane along the lambda ray

def resolution(des,c,EPD):
    cen=[np.mean(des.trace(l,c,npup=1200,Csys=CS)[0]) for l in BAND]; rows=[]
    for i,l in enumerate(BAND):
        U,_=des.trace(l,c,npup=2500,Csys=CS); rms=(U-U.mean()).std()*1e3
        i0,i1=max(0,i-1),min(len(BAND)-1,i+1); P=abs((lam_nm[i1]-lam_nm[i0])/(cen[i1]-cen[i0]))
        diff=1.028*l*(plane_r(l)/EPD)
        dl=np.hypot(np.hypot(SLIT,2.355*rms),np.hypot(PIX,diff))/1e3*P
        floor=np.hypot(SLIT,np.hypot(PIX,diff))/1e3*P
        rows.append(dict(lam=int(l*1e3),rms=float(rms),dl=float(dl),floor=float(floor)))
    return rows, abs(cen[-1]-cen[0])

SCALE=np.array([0.1,0.01,1e-4])
def band_cost(Rs,c,EPD):
    des=R.Design("P0400",N,ALPHA,rA,rA,Rs,Rs,EPD,BAND); out=[]
    for l in BAND:
        U,_=des.trace(l,list(c),npup=300,Csys=CS)
        if len(U)<50: return 1e6
        out.append((U-U.mean()).std()*1e3)
    out=np.array(out); return out.mean()+0.6*out.max()
def solve(EPD, Rs=None):
    """Fit c1,c2,c3 on the fixed datasheet plane (R fixed from the blank). If Rs is None, R is also
    freed -- used only as a cross-check on the blank-derived radius."""
    free_R = Rs is None
    def obj(z):
        if free_R: return band_cost(100+z[0], z[1:]*SCALE, EPD)
        return band_cost(Rs, z*SCALE, EPD)
    seeds=([1.8,-1.8,0],[-2,-1.8,0],[0,0,0])
    if free_R: seeds=[[0,*s_] for s_ in seeds]
    best=min((minimize(obj,s_,method="Nelder-Mead",options=dict(xatol=1e-4,fatol=1e-3,maxiter=2000)) for s_ in seeds),
             key=lambda r:r.fun)
    if free_R: return 100+best.x[0], list(best.x[1:]*SCALE)
    return Rs, list(best.x*SCALE)

def focal_distance(des,c,l):
    """distance along the lambda ray where that wavelength is sharpest (for the focal-curve plot)"""
    def s(d):
        dc=np.array([0,np.sin(beta(l)),np.cos(beta(l))]); ex=np.array([1.,0,0]); ed=np.cross(-dc,ex); ed/=np.linalg.norm(ed)
        U,_=des.trace(l,c,npup=300,Csys=(d*dc,-dc,ed,ex)); return (U-U.mean()).std()
    return minimize_scalar(s,bounds=(85,112),method="bounded").x

def main():
    print("="*78); print("  VALIDATION vs Shimadzu P0400-01 (spherical, flat-field) -- detector on the datasheet plane"); print("="*78)
    print(f"\n[1] beta at 400 nm {np.degrees(beta(0.4)):+.2f} deg, 1100 nm {np.degrees(beta(1.1)):+.2f} deg")
    print(f"[2] datasheet detector plane (r'i, beta_i): perpendicular distance {PDIST:.2f} mm, normal at beta {PBETA:+.2f} deg")
    print(f"    -> plane distance at 400 nm {plane_r(0.4):.2f} mm, 1100 nm {plane_r(1.1):.2f} mm (datasheet 97.1 / 99.5)")
    print(f"[3] radius of curvature from the blank sag: R = {RSUB:.2f} mm (fixed input)")
    Rfree,_=solve(rA/5.3,None)
    print(f"    cross-check: freeing R in the ray-trace fit gives R = {Rfree:.2f} mm")
    out={"R_blank":RSUB,"R_free":Rfree}
    for EPD,tag in ((rA/5.3,"f/5.3"),(rA/10,"f/10")):
        Rs,c=solve(EPD,RSUB); des=R.Design("P0400",N,ALPHA,rA,rA,Rs,Rs,EPD,BAND)
        rows,L=resolution(des,c,EPD)
        m=np.mean([r['dl'] for r in rows]); w=max(r['dl'] for r in rows)
        print(f"\n[{tag}] R = {Rs:.2f} mm (blank) | solved c1={c[0]:+.4f} c2={c[1]:+.5f} c3={c[2]:+.2e}")
        print(f"      detector length (traced) = {L:.2f} mm  | from dispersion {(1100-400)/PUB['disp']:.2f} | datasheet {PUB['L']}")
        print("      dlam  "+"  ".join(f"{r['lam']}:{r['dl']:.2f}" for r in rows)+f"   mean {m:.2f} worst {w:.2f}")
        print("      floor "+"  ".join(f"{r['floor']:.2f}" for r in rows))
        out[tag]=dict(R=Rs,c=c,L=L,rows=rows,mean=m,worst=w)
    # focal curve at f/5.3 solution vs the datasheet plane
    Rs,c=out["f/5.3"]["R"],out["f/5.3"]["c"]; des=R.Design("P0400",N,ALPHA,rA,rA,Rs,Rs,rA/5.3,BAND)
    fc=[focal_distance(des,c,l) for l in BAND]; out["focal"]=fc
    c=out["f/5.3"]["c"]; ypk=-c[0]/(2*c[1]); nu=R.nu_of(np.array([-20.,0.,ypk,20.]),N,*c)
    print(f"\n    groove density (f/5.3 solution): {nu[0]:.2f} at y=-20 | {nu[1]:.2f} at the POLE | "
          f"peak {nu[2]:.2f} at y={ypk:+.1f} mm | {nu[3]:.2f} at y=+20")
    out["groove"]=dict(y_peak=ypk,nu_minus20=nu[0],nu_pole=nu[1],nu_peak=nu[2],nu_plus20=nu[3])
    groove_figure(c)
    json.dump(out,open("validate_p0400.json","w"),indent=1,default=float)
    make_figure(out); print("\nsaved validate_p0400.png and validate_p0400.json")

def make_figure(out):
    CA,CB,CG="#1aa5cf","#e8912a","#2f8f5b"
    fig,(a1,a2,a3)=plt.subplots(1,3,figsize=(15.5,4.6))
    fig.suptitle("Validation vs Shimadzu P0400-01 (spherical flat-field) — detector fixed on the datasheet plane",
                 fontsize=12.5,fontweight="bold",y=1.02)
    L=out["f/5.3"]["L"]
    lams=[r['lam'] for r in out["f/5.3"]["rows"]]
    a1.axhline(PUB["L"],color=CG,ls="--",lw=1.2); a1.text(0.5,PUB["L"]+0.6,f"datasheet L = {PUB['L']} mm",color=CG,fontsize=9,ha="center")
    a1.bar([0],[L],width=0.5,color=CA); a1.bar([1],[(1100-400)/PUB['disp']],width=0.5,color=CB)
    a1.set_xticks([0,1]); a1.set_xticklabels([f"ray trace\n{L:.2f} mm",f"from dispersion\n{(1100-400)/PUB['disp']:.2f} mm"])
    a1.set_xlim(-0.6,1.6); a1.set_ylim(0,32); a1.set_ylabel("detector length (mm)"); a1.set_title("Detector length",fontsize=10.5)
    a2.plot(lam_nm,[plane_r(l) for l in BAND],"-",color=CG,lw=2.2,label="datasheet detector plane")
    a2.plot(lam_nm,out["focal"],"s",color=CB,ms=7,label="best focus, solved grooves")
    a2.plot([400,1100],[97.1,99.5],"o",mfc="none",mec=CG,ms=11,mew=1.8,label="datasheet r'₁ / r'₂")
    a2.set_xlabel("wavelength (nm)"); a2.set_ylabel("distance from grating pole (mm)")
    a2.set_title("Focus vs the flat detector plane\n(crossings = defocus nulls)",fontsize=10.5); a2.grid(alpha=0.3); a2.legend(fontsize=8)
    for tag,col,mk in (("f/5.3",CA,"o"),("f/10",CB,"s")):
        a3.plot(lams,[r['dl'] for r in out[tag]["rows"]],mk+"-",color=col,lw=2,label=f"{tag}  (mean {out[tag]['mean']:.2f} nm)")
    a3.plot(lams,[r['floor'] for r in out["f/10"]["rows"]],":",color="0.4",lw=1.4,label="slit + pixel + diffraction floor")
    a3.set_xlabel("wavelength (nm)"); a3.set_ylabel("resolution Δλ (nm)"); a3.set_ylim(0,3.3)
    a3.set_title("Resolution on the flat plane\n(floor reached at the two null wavelengths)",fontsize=10.5)
    a3.grid(alpha=0.3); a3.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("validate_p0400.png",dpi=150,bbox_inches="tight"); plt.close(fig)

def groove_figure(c):
    y=np.linspace(-20,20,401); nu=R.nu_of(y,N,*c); ypk=-c[0]/(2*c[1]); nupk=R.nu_of(np.array([ypk]),N,*c)[0]
    fig,ax=plt.subplots(figsize=(8.5,4.4))
    ax.plot(y,nu,color="#1aa5cf",lw=2.4)
    ax.axvline(0,color="0.6",lw=0.8,ls=":"); ax.axhline(N,color="0.75",lw=0.8,ls=":")
    ax.plot([0],[N],"o",color="#2f8f5b",ms=8); ax.annotate(f"pole: {N:.0f} l/mm (datasheet value)",(0,N),(-18.5,N+0.9),
                fontsize=9,color="#2f8f5b",arrowprops=dict(arrowstyle="->",color="#2f8f5b"))
    ax.plot([ypk],[nupk],"s",color="#e8912a",ms=7); ax.annotate(f"peak {nupk:.2f} l/mm at y = {ypk:+.1f} mm",(ypk,nupk),(ypk+1.2,nupk-3.2),
                fontsize=9,color="#e8912a",arrowprops=dict(arrowstyle="->",color="#e8912a"))
    ax.set_xlabel("position across the 40 mm blank, y (mm)   [y > 0 toward the entrance-slit side]")
    ax.set_ylabel("local groove density ν (l/mm)")
    ax.set_title(f"P0400-01 solved groove profile: ν(y) = 400 {c[0]:+.4f}y {c[1]:+.5f}y² {c[2]:+.1e}y³",fontsize=10.5)
    ax.grid(alpha=0.3); fig.tight_layout(); fig.savefig("p0400_groove_profile.png",dpi=150); plt.close(fig)

if __name__=="__main__":
    main()
