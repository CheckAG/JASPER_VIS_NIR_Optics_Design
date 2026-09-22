#!/usr/bin/env python3
"""Two published spot diagrams from the shared standalone VLS ray tracer (vls_raytrace.py):
   Panel A  f/5.3  SPHERICAL  — wide open: bright, aberration-limited (the tradeoff's fast end)
   Panel B  f/8    TOROIDAL   — stopped + astigmatism nulled: sharp AND ~full throughput

Only two knobs change between panels:
   EPD  -> f-number   (f/# = rA/EPD)     : the high->medium RESOLUTION knob
   rho  -> substrate  (Ry=spherical; rho<Ry=toroid) : the THROUGHPUT knob (astig line height)
Everything else (tracer, grating eq, best-focus) is identical and untouched.
"""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import vls_raytrace as R

SLIT=20.0; PIX=8.0
band=[0.5,0.6,0.7,0.8,0.9,1.0]
lam_nm=np.array(band)*1e3
COL={500:"#1aa5cf",600:"#e8912a",700:"#e34b3a",800:"#b23a44",900:"#8f4f63",1000:"#6f6480"}

# ---- Simulated custom VLS grating (Phase 1): 230 l/mm, alpha 12.34, arms 137.357/126.75, Ry 133 ----
N=230.0; ALPHA=12.34; rA=137.357; rB=126.75; Ry=133.0
# toroid sag radius that nulls astigmatism at band centre: rho = (cos a + cos b)/(1/rA + 1/rB)
a=np.radians(ALPHA)
bc=np.arcsin((0.75*1000)/(1e6/N)-np.sin(a))
RHO_TOROID=(np.cos(a)+np.cos(bc))/(1.0/rA+1.0/rB)
print(f"toroid sag radius rho = {RHO_TOROID:.2f} mm  (spherical Ry = {Ry} mm)")

def fnum_to_EPD(fn): return rA/fn

PANELS=[
    dict(tag="A", fn=5.3, rho=Ry,        sub="SPHERICAL", note="wide open — bright, aberration-limited"),
    dict(tag="B", fn=8.0, rho=RHO_TOROID,sub="TOROIDAL",  note="stopped + astigmatism nulled — sharp & ~full throughput"),
]

def run(fn, rho):
    EPD=fnum_to_EPD(fn)
    d=R.Design("sim",N,ALPHA,rA,rB,Ry,rho,EPD,band)
    c=[0,0,0]                                    # classical grooves = no unknowns
    rBf=d.best_focus(c); Cs=d.detector(c,rBf)
    cen=[np.mean(d.trace(l,c,npup=1500,Csys=Cs)[0]) for l in band]
    rows=[]
    for i,l in enumerate(band):
        U,V=d.trace(l,c,npup=2500,Csys=Cs); cu=U.mean()
        aa=max(0,i-1); bb=min(5,i+1)
        plate=abs((lam_nm[bb]-lam_nm[aa])/(cen[bb]-cen[aa]))
        rms=(U-cu).std()*1000; line=(V.max()-V.min())
        diff=1.028*l*(d.rB/EPD)
        dl=np.hypot(np.hypot(SLIT,2.355*rms),np.hypot(PIX,diff))/1000*plate
        rows.append(dict(l=l,U=U-cu,V=V,rms=rms,line=line,plate=plate,dl=dl))
    return EPD,rBf,rows

fig,axs=plt.subplots(2,6,figsize=(15,7.0))
for prow,P in enumerate(PANELS):
    EPD,rBf,rows=run(P["fn"],P["rho"])
    meandl=np.mean([r["dl"] for r in rows]); meanline=np.mean([r["line"] for r in rows])
    print(f"\nPanel {P['tag']}: f/{P['fn']} {P['sub']}  EPD={EPD:.1f}mm r'={rBf:.1f}mm  "
          f"mean dlam {meandl:.2f} nm  astig line {meanline:.2f} mm")
    for i,r in enumerate(rows):
        ax=axs[prow,i]; k=int(r["l"]*1000)
        ax.scatter(r["U"]*1000,r["V"],s=1.6,color=COL[k],alpha=0.5,lw=0)
        ax.axvspan(-SLIT/2,SLIT/2,color=COL[k],alpha=0.10)
        ax.set_xlim(-140,140); ax.set_ylim(-1.6,1.6)
        ax.set_title(f"{k} nm\nΔλ={r['dl']:.2f} nm",fontsize=9.5,color=COL[k])
        ax.axhline(0,color="0.85",lw=.6); ax.axvline(0,color="0.85",lw=.6); ax.tick_params(labelsize=7)
        if i==0: ax.set_ylabel(f"Panel {P['tag']}  f/{P['fn']}  {P['sub']}\ncross-dispersion (mm)",fontsize=8.5)
        if prow==1: ax.set_xlabel("dispersion (µm)",fontsize=8)
        ax.text(0,-1.45,f"RMS {r['rms']:.0f}µm · line {r['line']:.2f}mm",fontsize=6.5,ha="center",color="0.35")
    axs[prow,5].text(1.02,0.5,f"mean Δλ {meandl:.2f} nm\nastig line {meanline:.2f} mm",
                     transform=axs[prow,5].transAxes,fontsize=8.5,va="center",color="0.25")

fig.suptitle("JASPER simulated VLS grating (230 l/mm) — ray-traced spots: the resolution/throughput tradeoff\n"
             "Panel A f/5.3 spherical (bright)  vs  Panel B f/8 toroidal (sharp + full throughput) · standalone tracer, classical grooves, best focus",
             fontsize=10.5,y=1.02)
fig.text(0.5,-0.02,"Shaded band = 20 µm slit. Between panels only EPD (f-number, via the aperture stop) and rho (spherical Ry → toroid) change; "
         "the tracer, grating equation and best-focus search are identical. Line = astigmatic height (throughput proxy).",
         fontsize=8,color="0.4",ha="center")
fig.tight_layout()
fig.savefig("jasper_two_panel_spots.png",dpi=150,bbox_inches="tight")
print("\nwrote jasper_two_panel_spots.png")
