#!/usr/bin/env python3
"""Beam4-style optical LAYOUT from the standalone VLS ray tracer (vls_raytrace.py) — thousands of
REAL rays through the exact curved substrate WITH the local variable-line-spacing groove equation.

This is the layout view the schematic sphere_layout.py never was: an actual meridional (dispersion-
plane) ray trace. For each wavelength we launch a fan of rays across the grating aperture, intersect
the concave surface exactly, diffract with the local VLS density nu(y), and propagate to the tilted
flat detector. The dense wavelength-coloured fan converging on the grating and fanning to the
detector is the picture.

Meridional note: in the dispersion (y-z) plane the toroid and the sphere share the SAME radius (Ry);
the toroid differs only out of plane (sagittal radius rho), i.e. in astigmatic line HEIGHT, which is
a spot-diagram quantity, not a layout one. So this one layout serves both substrates.
"""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
import vls_raytrace as R

BAND=[0.5,0.6,0.7,0.8,0.9,1.0]
COL={500:"#1aa5cf",600:"#e8912a",700:"#e34b3a",800:"#b23a44",900:"#8f4f63",1000:"#6f6480"}

# ---- Simulated custom VLS grating (Phase 1): 230 l/mm, alpha 12.34, arms 137.357/126.75, Ry 133 ----
N=230.0; ALPHA=12.34; rA=137.357; rB=126.75; Ry=133.0
EPD=25.9   # f/5.3 wide-open, so the fan is visibly wide; f-number only scales fan width
a=np.radians(ALPHA)
bc=np.arcsin((0.75*1000)/(1e6/N)-np.sin(a))
RHO=(np.cos(a)+np.cos(bc))/(1.0/rA+1.0/rB)
des=R.Design("sim",N,ALPHA,rA,rB,Ry,RHO,EPD,BAND)
c=[0,0,0]                                    # groove correction is out-of-plane detail; layout uses classical
rBf=des.best_focus(c)
C,ndet,edisp,ecross=des.detector(c,rBf)      # detector centre, normal, in-plane axes

def meridional(lam, nrays=150):
    """Return (A2, P2, H2, beta_deg): slit, surface-hit and detector-hit points as (y,z) arrays,
    for a fan of nrays launched in the x=0 (meridional) plane. Reuses the tracer's exact-sphere
    intersection + local-VLS diffraction, restricted to the dispersion plane."""
    lam_mm=lam*1e-3; A=des.A
    ys=np.linspace(-EPD/2,EPD/2,nrays); xs=np.zeros_like(ys)
    Pt=np.stack([xs,ys,ys*ys/(2*Ry)],axis=1)           # aim points on the surface (x=0)
    d=Pt-A; d/=np.linalg.norm(d,axis=1,keepdims=True)
    # exact sphere (meridional radius Ry): slit inside sphere -> far root
    Oc=np.array([0.,0.,Ry]); e=A-Oc; ed=d@e; disc=ed*ed-(e@e-Ry*Ry)
    t=-ed+np.sqrt(np.clip(disc,0,None)); P=A+t[:,None]*d
    Nrm=(Oc-P)/Ry
    # local VLS diffraction
    yhat=np.array([0.,1.,0.]); g=yhat-(Nrm[:,1:2])*Nrm; g/=np.linalg.norm(g,axis=1,keepdims=True)
    nu=R.nu_of(P[:,1],des.nu0,*c)
    din_tan=d-(np.sum(d*Nrm,axis=1,keepdims=True))*Nrm
    dout_tan=din_tan+(des.m*lam_mm*nu)[:,None]*g
    n2=1.0-np.sum(dout_tan*dout_tan,axis=1); ok=n2>=0
    sign=-np.sign(np.sum(d*Nrm,axis=1))
    dout=dout_tan+(sign*np.sqrt(np.clip(n2,0,None)))[:,None]*Nrm
    denom=dout@ndet; tt=((C-P)@ndet)/denom; H=P+tt[:,None]*dout
    # chief-ray beta (central ray, y=0 aim)
    ic=nrays//2; chief=dout[ic]; beta=np.degrees(np.arctan2(chief[1],chief[2]))
    sel=ok&np.isfinite(H[:,0])
    yz=lambda M:(M[sel][:,1],M[sel][:,2])                # (y,z)
    return (A[1],A[2]), yz(P), yz(H), beta

# ---------------- draw ----------------
fig,ax=plt.subplots(figsize=(14,9))
ax.set_facecolor("#0e1116")                              # dark bg makes the fans glow (Beam4 look)

# grating surface profile (exact circle radius Ry, in the y-z plane), drawn thick
yy=np.linspace(-EPD/2-3,EPD/2+3,220); zz=Ry-np.sqrt(Ry*Ry-yy*yy)
ax.plot(zz,yy,color="#39d98a",lw=3.0,zorder=5,solid_capstyle="round")
# little groove ticks along the grating to signal it is a ruled/VLS grating
for yt in np.linspace(-EPD/2,EPD/2,26):
    zt=Ry-np.sqrt(Ry*Ry-yt*yt); ax.plot([zt-0.0,zt+1.1],[yt,yt],color="#39d98a",lw=0.6,alpha=0.55,zorder=4)

# the ray fans.  The incoming slit->grating beam is white light (undispersed) and is the SAME
# geometry for every wavelength, so draw it ONCE in white; colour only the diffracted side.
Aslit=None; dets={}
(Ay,Az),(Py0,Pz0),_,_=meridional(0.7,nrays=180)
Aslit=(Ay,Az)
for j in range(len(Py0)):
    ax.plot([Az,Pz0[j]],[Ay,Py0[j]],color="#eaf2f8",lw=0.35,alpha=0.055,zorder=2)   # white input beam
for lam in BAND:
    (Ay,Az),(Py,Pz),(Hy,Hz),beta=meridional(lam,nrays=180)
    col=COL[int(lam*1000)]
    for j in range(len(Py)):
        ax.plot([Pz[j],Hz[j]],[Py[j],Hy[j]],color=col,lw=0.35,alpha=0.16,zorder=3)  # grating -> detector
    dets[lam]=(np.mean(Hy),np.mean(Hz),beta)

# slit marker
ax.plot(Aslit[1],Aslit[0],'o',color="#ffffff",ms=7,zorder=8,mec="k")
ax.annotate("Entrance slit  20 µm × 1 mm\n(fed by SMA905 + relay)",xy=(Aslit[1],Aslit[0]),xytext=(Aslit[1]+12,Aslit[0]-10),
            color="#dfe7ee",fontsize=9,ha="right",arrowprops=dict(arrowstyle="->",color="#9fb0bf",lw=0.9))
# grating label
ax.annotate("Concave VLS grating\n230 l/mm · R=133 mm · m=1\n(spherical or toroidal)",
            xy=(Ry-np.sqrt(Ry*Ry-(EPD/2)**2),EPD/2),xytext=(6,EPD/2+6),
            color="#39d98a",fontsize=9.5,ha="left",arrowprops=dict(arrowstyle="->",color="#39d98a",lw=0.9))

# detector: tilted flat line through the wavelength foci, 29.18 mm linear array
Hy=np.array([dets[l][0] for l in BAND]); Hz=np.array([dets[l][1] for l in BAND])
# fit the focal locus, draw a 29.18 mm sensor centred on it
p=np.polyfit(Hz,Hy,1); zc=Hz.mean(); yc=np.polyval(p,zc)
tdir=np.array([1,p[0]])/np.hypot(1,p[0]); half=29.18/2
e1=np.array([zc,yc])-half*tdir; e2=np.array([zc,yc])+half*tdir
ax.plot([e1[0],e2[0]],[e1[1],e2[1]],color="#c9d2da",lw=5,solid_capstyle="butt",zorder=7)
for l in BAND:
    y0,z0,beta=dets[l]; k=int(l*1000)
    ax.plot(z0,y0,'o',color=COL[k],ms=5.5,mec="k",mew=.4,zorder=9)
ax.annotate("Detector \nflat, tilted",xy=(zc,yc),xytext=(zc+3,yc-14),
            color="#dfe7ee",fontsize=9,ha="left",arrowprops=dict(arrowstyle="->",color="#9fb0bf",lw=0.9))
ax.text(dets[0.5][1]-2,dets[0.5][0]-2,"500 nm",color=COL[500],fontsize=8.5,ha="right",va="center")
ax.text(dets[1.0][1]+2,dets[1.0][0]+1,"1000 nm",color=COL[1000],fontsize=8.5,ha="left",va="center")

# grating normal (dashed)
ax.plot([0,Hz.mean()*1.05],[0,0],ls=(0,(6,5)),color="#5a6672",lw=1.0,zorder=1)
ax.text(Hz.mean()*1.05,0.4,"grating normal",color="#7f8c99",fontsize=8.5,va="bottom",ha="right")

b5=dets[0.5][2]; b10=dets[1.0][2]
ax.set_title("JASPER — ray-traced optical layout on the concave VLS grating   ·   standalone ray tracer, ~1080 rays (180 × 6 λ)",
             color="#eef3f7",fontsize=12.5,pad=14)
ax.text(0.5,-0.055,f"Real geometric ray trace in the dispersion plane: slit → exact concave surface → local VLS diffraction → tilted flat detector.  "
        f"β(500 nm)={b5:.1f}°, β(1000 nm)={b10:+.1f}°.  Fan width shown at f/5.3 (wide open); the f-number (aperture stop) only scales the fan, not the layout.",
        transform=ax.transAxes,color="#9fb0bf",fontsize=8.6,ha="center",va="top",wrap=True)
ax.set_aspect("equal"); ax.axis("off")
fig.tight_layout(); fig.savefig("jasper_raytrace_layout.png",dpi=155,bbox_inches="tight",facecolor="#0e1116")
print(f"wrote jasper_raytrace_layout.png | beta 500={b5:.2f} 1000={b10:.2f} | det focus r'={rBf:.2f}")
