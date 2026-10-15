#!/usr/bin/env python3
"""Per-wavelength spot diagrams from the light-path-function model, for the JASPER designs.
IDEALISED BOUND (F20 nulled at every wavelength -- as if the detector followed the focal curve; a single
groove chirp on a flat detector cannot do this): dispersion spread = coma (F30) + 20um slit; cross-dispersion
height = astigmatism (F02). Same theory as the resolution numbers -> consistent with them, no
fragile ray trace. Spherical blank shown (toroidal collapses the height)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
SLIT=20.0; PIX=8.0
COL={500:"#1aa5cf",600:"#e8912a",700:"#e34b3a",800:"#b23a44",900:"#8f4f63",1000:"#6f6480"}
WAVES=[0.5,0.6,0.7,0.8,0.9,1.0]

def make(name, tag, N, alpha_deg, rA, rB, Rtan, rho, EPD, bref_deg, outfile, ylim, astig_txt):
    d_nm=1e6/N; a=np.radians(alpha_deg); bi=np.radians(bref_deg); m=1
    def beta(l): return np.arcsin(m*(l*1000)/d_nm-np.sin(a))
    def rp(b): return rB/np.cos(b-bi)
    def plate(b,r): return d_nm*np.cos(b)*np.cos(b-bi)/(m*r)   # Shimadzu eq (24)
    def F30(b,r): return (np.sin(a)/rA)*(np.cos(a)**2/rA-np.cos(a)/Rtan)+(np.sin(b)/r)*(np.cos(b)**2/r-np.cos(b)/Rtan)
    def F02(b,r): return 1.0/rA+1.0/r-(np.cos(a)+np.cos(b))/rho
    rng=np.random.default_rng(3)
    th=rng.uniform(0,2*np.pi,1400); rr=np.sqrt(rng.uniform(0,1,1400))*EPD/2
    W=rr*np.cos(th); L=rr*np.sin(th); slit=rng.uniform(-SLIT/2,SLIT/2,1400)
    fig,axs=plt.subplots(1,6,figsize=(15,3.5))
    for ax,w in zip(axs,WAVES):
        b=beta(w); r=rp(b); P=plate(b,r)
        dy=(r/np.cos(b))*1.5*F30(b,r)*W**2*1000 + slit; dy-=dy.mean()
        dx=r*F02(b,r)*L*1000
        rmsc=np.sqrt(np.mean(((r/np.cos(b))*1.5*F30(b,r)*W**2*1000)**2))
        dl=np.hypot(np.hypot(SLIT,2.355*rmsc),PIX)/1000*P
        ax.scatter(dy,dx/1000.0,s=2.2,color=COL[int(w*1000)],alpha=0.55,lw=0)
        ax.axvspan(-SLIT/2,SLIT/2,color=COL[int(w*1000)],alpha=0.10)
        ax.set_xlim(-70,70); ax.set_ylim(-ylim,ylim)
        ax.set_title(f"{int(w*1000)} nm\nΔλ={dl:.2f} nm",fontsize=9.5,color=COL[int(w*1000)])
        ax.axhline(0,color="0.8",lw=0.6); ax.axvline(0,color="0.8",lw=0.6); ax.tick_params(labelsize=7)
        if ax is axs[0]: ax.set_ylabel("cross-dispersion (mm)",fontsize=8.5)
        ax.set_xlabel("dispersion (µm)",fontsize=8)
        ax.text(0.5,-ylim*0.9,f"RMS$_x$ {rmsc:.0f}µm",fontsize=7,ha="center",color="0.35")
    fig.suptitle(f"{name} — spot diagrams (light-path-function, IDEALISED BOUND (defocus nulled at every λ — not reachable on a flat detector), f/5.3, spherical blank)  ·  {tag}",fontsize=10.5,y=1.06)
    fig.text(0.5,-0.02,f"Shaded band = 20 µm slit image. Dispersion axis (horizontal) sets resolution; cross-dispersion height (vertical) is astigmatism = throughput. {astig_txt}",fontsize=8,color="0.4",ha="center")
    fig.tight_layout(); fig.savefig(outfile,dpi=155,bbox_inches="tight"); print("wrote",outfile)

if __name__=="__main__":
    bc=np.arcsin((0.75*1000)/(1e6/230.0)-np.sin(np.radians(12.34)))
    make("SIMULATED DESIGN (concave VLS, spherical)","each wavelength: sharp in dispersion (~0.9 nm), ~1.0 mm astigmatic line",
         230,12.34,137.357,126.75,133.0,133.0,26.0,-2.4,"jasper_simulated_spots.png",1.4,
         "A toroidal blank collapses the height to ~0.04 mm.")
    make("DESIGN C — Shimadzu P0400-01","each wavelength: sharp in dispersion (~0.6 nm), ~2.3 mm astigmatic line",
         400,20.0,96.9,96.93,99.61,99.61,18.3,-7.51,"jasper_p0400_spots.png",1.4,
         "A toroidal blank (Shimadzu's toroidal family) collapses the height to ~0.15 mm.")
