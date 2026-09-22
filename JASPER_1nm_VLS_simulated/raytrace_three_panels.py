#!/usr/bin/env python3
"""Three published spot diagrams from the shared standalone VLS ray tracer (vls_raytrace.py):
   Panel A  f/5.3  SPHERICAL  classical   — wide open: bright, aberration-limited (fast end)
   Panel B  f/8    TOROIDAL   classical   — stopped + astig nulled: sharper & ~full throughput
   Panel C  f/10   TOROIDAL   VLS-corrected — the groove correction kills the blue-edge aberration -> ~1 nm

Knobs changed between panels (nothing in the tracer changes):
   EPD  -> f-number (f/# = rA/EPD)                 : shrinks aberration RMS (high->medium res)
   rho  -> substrate (Ry=spherical; rho<Ry=toroid) : throughput (astig line height), ~0 res effect
   c    -> VLS groove profile ([0,0,0]=classical; fitted=corrected) : nulls FIELD aberration at edges
The A->B->C path is exactly "what it takes to reach true ~1 nm across 500-1000 nm."
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
a=np.radians(ALPHA)
bc=np.arcsin((0.75*1000)/(1e6/N)-np.sin(a))
RHO_TOROID=(np.cos(a)+np.cos(bc))/(1.0/rA+1.0/rB)
print(f"toroid sag radius rho = {RHO_TOROID:.2f} mm  (spherical Ry = {Ry} mm)")

def fnum_to_EPD(fn): return rA/fn

PANELS=[
    dict(tag="A", fn=5.3, rho=Ry,        corr=False, sub="spherical, classical",
         note="wide open — bright, aberration-limited"),
    dict(tag="B", fn=8.0, rho=RHO_TOROID,corr=False, sub="toroidal, classical",
         note="stopped + astig nulled — sharper & ~full throughput"),
    dict(tag="C", fn=10.0,rho=RHO_TOROID,corr=True,  sub="toroidal, VLS-corrected",
         note="groove correction clears the blue edge — ~1 nm across the band (best case)"),
]

def run(fn, rho, corr):
    EPD=fnum_to_EPD(fn)
    d=R.Design("sim",N,ALPHA,rA,rB,Ry,rho,EPD,band)
    c=R.fit_vls(d) if corr else [0,0,0]
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
    return EPD,rBf,c,rows

fig,axs=plt.subplots(3,6,figsize=(15,10.2))
summary=[]
for prow,P in enumerate(PANELS):
    EPD,rBf,c,rows=run(P["fn"],P["rho"],P["corr"])
    meandl=np.mean([r["dl"] for r in rows]); worst=max(r["dl"] for r in rows)
    meanline=np.mean([r["line"] for r in rows])
    summary.append((P,EPD,rBf,c,meandl,worst,meanline,rows))
    print(f"\nPanel {P['tag']}: f/{P['fn']} {P['sub']}  EPD={EPD:.1f}mm r'={rBf:.1f}mm  "
          f"mean {meandl:.2f} worst {worst:.2f} nm  astig line {meanline:.2f} mm")
    for r in rows: print(f"    {int(r['l']*1e3):>4}nm  RMS {r['rms']:>4.0f}um  dlam {r['dl']:.2f} nm")
    for i,r in enumerate(rows):
        ax=axs[prow,i]; k=int(r["l"]*1000)
        ax.scatter(r["U"]*1000,r["V"],s=1.6,color=COL[k],alpha=0.5,lw=0)
        ax.axvspan(-SLIT/2,SLIT/2,color=COL[k],alpha=0.10)
        ax.set_xlim(-140,140); ax.set_ylim(-1.6,1.6)
        edge = " ←edge" if k==500 else ""
        ax.set_title(f"{k} nm\nΔλ={r['dl']:.2f} nm{edge}",fontsize=9.5,color=COL[k])
        ax.axhline(0,color="0.85",lw=.6); ax.axvline(0,color="0.85",lw=.6); ax.tick_params(labelsize=7)
        if i==0: ax.set_ylabel(f"Panel {P['tag']}  f/{P['fn']}\n{P['sub']}\ncross-disp (mm)",fontsize=8.2)
        if prow==2: ax.set_xlabel("dispersion (µm)",fontsize=8)
        ax.text(0,-1.45,f"RMS {r['rms']:.0f}µm · line {r['line']:.2f}mm",fontsize=6.3,ha="center",color="0.35")
    axs[prow,5].text(1.03,0.5,f"mean {meandl:.2f} nm\nworst {worst:.2f} nm\nastig {meanline:.2f} mm",
                     transform=axs[prow,5].transAxes,fontsize=8.2,va="center",color="0.25")

fig.suptitle("JASPER simulated VLS grating (230 l/mm) — the path to ~1 nm across 500–1000 nm\n"
             "A: f/5.3 spherical (bright)  →  B: f/8 toroidal (throughput)  →  C: f/10 toroidal + VLS correction (clears the 500 nm blue edge)\n"
             "standalone ray tracer · best focus · Panel C VLS = best case (assumes fabricated groove profile matches the fit)",
             fontsize=10.3,y=1.03)
fig.text(0.5,-0.015,"Shaded band = 20 µm slit. The binding constraint is the 500 nm blue edge: aberration-limited, not slit/plate-limited. "
         "Stopping down (A→B) and the toroid help, but only the VLS groove correction (B→C) pulls the blue edge under ~1 nm.",
         fontsize=8,color="0.4",ha="center")
fig.tight_layout()
fig.savefig("jasper_three_panel_spots.png",dpi=150,bbox_inches="tight")
print("\nwrote jasper_three_panel_spots.png")

# emit a compact JSON of the numbers so the design record can cite them exactly
import json
out={P["tag"]:dict(fn=P["fn"],sub=P["sub"],EPD=round(EPD,1),rp=round(rBf,1),
                   c=[float(x) for x in (c if isinstance(c,(list,np.ndarray)) else [0,0,0])],
                   mean=round(md,2),worst=round(wo,2),astig_mm=round(ln,2),
                   per=[dict(lam=int(r["l"]*1e3),dl=round(r["dl"],2),rms=round(r["rms"],0),line=round(r["line"],2)) for r in rr])
     for (P,EPD,rBf,c,md,wo,ln,rr) in summary}
json.dump(out, open("jasper_three_panel_numbers.json","w"), indent=1)
print("wrote jasper_three_panel_numbers.json")
