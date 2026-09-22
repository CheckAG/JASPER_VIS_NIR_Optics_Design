#!/usr/bin/env python3
"""Light-background traced-ray layout in the style of the old vls_fa_figs.py 'Design A' figure,
but driven by the CURRENT standalone tracer (vls_raytrace.py) and the CURRENT consolidated design
(230 l/mm, R=133 mm, rA=137.357, rB=126.75, alpha=12.34). Real traced rays: a red
diverging entrance cone from the slit to the grating aperture, then wavelength-coloured diffracted
rays to a tilted flat detector placed on the true focal curve. Honest f-number-tradeoff numbers in
the info box (not the old single best-case figure)."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import vls_raytrace as R

BAND=[0.5,0.6,0.7,0.8,0.9,1.0]
COL={500:"#1aa5cf",600:"#e8912a",700:"#e34b3a",800:"#b23a44",900:"#8f4f63",1000:"#6f6480"}; RED="#d81f1f"
N=230.0; ALPHA=12.34; rA=137.357; rB=126.75; Ry=133.0; EPD=25.9
a=np.radians(ALPHA)
bc=np.arcsin((0.75*1000)/(1e6/N)-np.sin(a))
RHO=(np.cos(a)+np.cos(bc))/(1.0/rA+1.0/rB)
des=R.Design("sim",N,ALPHA,rA,rB,Ry,RHO,EPD,BAND)
c=[0,0,0]; rBf=des.best_focus(c); C,ndet,edisp,ecross=des.detector(c,rBf)

def meridional(lam,nrays=15):
    lam_mm=lam*1e-3; A=des.A
    ys=np.linspace(-EPD/2,EPD/2,nrays); xs=np.zeros_like(ys)
    Pt=np.stack([xs,ys,ys*ys/(2*Ry)],axis=1); d=Pt-A; d/=np.linalg.norm(d,axis=1,keepdims=True)
    Oc=np.array([0.,0.,Ry]); e=A-Oc; ed=d@e; disc=ed*ed-(e@e-Ry*Ry)
    t=-ed+np.sqrt(np.clip(disc,0,None)); P=A+t[:,None]*d; Nrm=(Oc-P)/Ry
    yhat=np.array([0.,1.,0.]); g=yhat-(Nrm[:,1:2])*Nrm; g/=np.linalg.norm(g,axis=1,keepdims=True)
    nu=R.nu_of(P[:,1],des.nu0,*c); din=d-(np.sum(d*Nrm,axis=1,keepdims=True))*Nrm
    dt=din+(des.m*lam_mm*nu)[:,None]*g; n2=1.0-np.sum(dt*dt,axis=1); ok=n2>=0
    sign=-np.sign(np.sum(d*Nrm,axis=1)); dout=dt+(sign*np.sqrt(np.clip(n2,0,None)))[:,None]*Nrm
    denom=dout@ndet; tt=((C-P)@ndet)/denom; H=P+tt[:,None]*dout
    ic=nrays//2; beta=np.degrees(np.arctan2(dout[ic,1],dout[ic,2]))
    sel=ok&np.isfinite(H[:,0]); yz=lambda M:np.stack([M[sel][:,2],M[sel][:,1]],axis=1)
    Dyz=yz(dout); Dyz/=np.linalg.norm(Dyz,axis=1,keepdims=True)
    return np.array([des.A[2],des.A[1]]), yz(P), yz(H), beta, Dyz

def focus_point(Pyz,Dyz):
    M=np.zeros((2,2)); b=np.zeros(2)
    for o,u in zip(Pyz,Dyz): Pn=np.eye(2)-np.outer(u,u); M+=Pn; b+=Pn@o
    return np.linalg.solve(M,b)

# gather geometry
Aslit,Pin,_,_,_=meridional(0.7,nrays=15)
foci=[]; betas={}
for l in BAND:
    _,P,H,beta,D=meridional(l,nrays=15); foci.append(focus_point(P,D)); betas[l]=beta
foci=np.array(foci)
p=np.polyfit(foci[:,0],foci[:,1],1); zc=foci[:,0].mean(); yc=np.polyval(p,zc)
tdir=np.array([1,p[0]])/np.hypot(1,p[0]); det_from_vert=abs(np.degrees(np.arctan2(tdir[0],tdir[1])))
plate=abs((1000-500)/np.hypot(*(foci[-1]-foci[0])))  # nm/mm along the sensor
spectrum_len=np.hypot(*(foci[-1]-foci[0])); dev=ALPHA+abs(0.5*(betas[0.5]+betas[1.0]))
fnum=rA/EPD

# ---------------- draw (light background, old style) ----------------
fig,ax=plt.subplots(figsize=(13.5,8.2))
xf=lambda y: y*y/(2*Ry)                                  # grating sag z(y)
AP=EPD/2
# grating dish (light blue) + groove ticks
yy=np.linspace(-AP-1.5,AP+1.5,80)
ax.fill(np.concatenate([xf(yy),[-6,-6]]),np.concatenate([yy,[AP+1.5,-AP-1.5]]),color="#cfe6f2",ec="k",lw=1,zorder=2)
ax.plot(xf(yy),yy,color="k",lw=2,zorder=3)
for yt in np.linspace(-AP,AP,15): ax.plot([xf(yt)-0.4,xf(yt)+0.4],[yt,yt],color="#3a6ea5",lw=0.5,zorder=3)
ax.text(-10,-AP+1,"Concave VLS grating\n(spherical / toroidal)",fontsize=9,ha="right",va="top")
# normal
ax.plot([0,44],[0,0],ls=(0,(7,5)),color="0.5",lw=1); ax.text(45,0,"normal",fontsize=9,va="center")
# entrance diverging cone (red) slit -> grating aperture
for pt in Pin: ax.plot([Aslit[0],pt[0]],[Aslit[1],pt[1]],color=RED,lw=0.6,alpha=0.7,zorder=1)
# slit (tilted, ⊥ entrance arm)
armdir=(-Aslit)/np.hypot(*Aslit); sperp=np.array([-armdir[1],armdir[0]])
s1=Aslit+2.0*sperp; s2=Aslit-2.0*sperp
ax.plot([s1[0],s2[0]],[s1[1],s2[1]],color="k",lw=5,solid_capstyle="butt",zorder=6)
ax.text(Aslit[0]+3,Aslit[1]+2,f"Entrance slit  (20 µm × 1 mm)\ntilted {ALPHA:.0f}° (⊥ entrance arm)",fontsize=9,ha="left")
ax.text(0.5*Aslit[0]+4,0.5*Aslit[1]+4,f"r$_A$={rA:.0f} mm",fontsize=9.5,style="italic",color=RED,
        rotation=np.degrees(np.arctan2(Aslit[1],Aslit[0])))
ax.add_patch(Arc((0,0),26,26,theta1=0,theta2=ALPHA,lw=1)); ax.text(15,2,rf"$\alpha$={ALPHA:.0f}°",fontsize=10)
# diffracted traced rays -> detector (colour per wavelength)
for l in BAND:
    _,P,H,_,_=meridional(l,nrays=15); col=COL[int(l*1000)]
    for j in range(len(P)): ax.plot([P[j,0],H[j,0]],[P[j,1],H[j,1]],color=col,lw=0.5,alpha=0.6,zorder=1)
# detector (flat, tilted, on the true focal curve)
half=spectrum_len/2+3; e1=np.array([zc,yc])-half*tdir; e2=np.array([zc,yc])+half*tdir
ax.plot([e1[0],e2[0]],[e1[1],e2[1]],color="k",lw=4,solid_capstyle="butt",zorder=6)
for l in BAND:
    f=foci[BAND.index(l)]; ax.plot(f[0],f[1],'o',color=COL[int(l*1000)],ms=5,mec="k",mew=0.3,zorder=7)
ax.text(foci[0][0]+3,foci[0][1]-1,"500 nm",color=COL[500],fontsize=9)
ax.text(foci[-1][0]+3,foci[-1][1],"1000 nm",color=COL[1000],fontsize=9,va="center")
ax.text(zc-2,foci[0][1]-6,f"Detector (flat, tilt {det_from_vert:.0f}°)\nr$_B$={rB:.0f} mm, L≈{spectrum_len:.0f} mm",
        fontsize=9,ha="center",va="top")
# info box (honest numbers)
box=(f"Finite-arm mount (real diverging slit), traced rays\n"
     f"concave VLS grating: R={Ry:.0f} mm, {N:.0f} l/mm, m=1  ·  spherical or toroidal substrate\n"
     f"α={ALPHA:.1f}°, deviation {dev:.0f}° (slit & detector separated)\n"
     f"r$_A$={rA:.0f}, r$_B$={rB:.0f} mm, f/{fnum:.1f}, plate {plate:.1f} nm/mm, 500–1000 nm → ~{spectrum_len:.0f} mm on a linear CCD\n"
     f"Δλ is an f-number tradeoff: ~2.3 nm at f/5.3 → ~1.4 nm at f/8 → ~1 nm at f/10 (VLS best case); floor ~0.6 nm\n"
     f"Throughput: spherical → ~1 mm astig line (~19% on pixel); toroidal → ~0.04 mm (~100%)")
# centred horizontally (x = axes fraction 0.5), vertical position in data coords
ax.text(0.5,Aslit[1]+26,box,transform=ax.get_yaxis_transform(),fontsize=8.2,family="monospace",
        va="top",ha="center",ma="center",bbox=dict(boxstyle="round,pad=0.5",fc="#f7f7f7",ec="#888"))
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("JASPER ~1 nm concave VLS grating — finite-arm layout, traced rays",fontsize=13)
ax.set_xlim(-48,rB+34); ax.set_ylim(foci[0][1]-16,Aslit[1]+34)
fig.tight_layout(); fig.savefig("jasper_layout_light.png",dpi=150,bbox_inches="tight",facecolor="white")
print(f"wrote jasper_layout_light.png | beta {betas[0.5]:.2f}->{betas[1.0]:.2f} | det tilt {det_from_vert:.1f} | "
      f"plate {plate:.1f} | spectrum {spectrum_len:.2f} | dev {dev:.1f}")
