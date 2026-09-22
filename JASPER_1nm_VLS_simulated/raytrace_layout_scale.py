#!/usr/bin/env python3
"""To-SCALE engineering layout of the JASPER concave-VLS spectrometer core, from the standalone
ray tracer (vls_raytrace.py). Everything is in true millimetres with equal aspect: the entrance arm
L_A, the diffracted arm L_H, the concave grating radius, the slit orientation (perpendicular to the
entrance arm) and the detector drawn on the ACTUAL focal locus at its true tilt. Dimension callouts
(L_A, L_H, alpha, beta span, detector tilt, spectrum length) are measured off the trace, not drawn
by hand. The faint ray fan is the same real VLS trace, dialled back so the geometry reads first.
"""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import vls_raytrace as R

BAND=[0.5,0.6,0.7,0.8,0.9,1.0]
COL={500:"#1aa5cf",600:"#e8912a",700:"#e34b3a",800:"#b23a44",900:"#8f4f63",1000:"#6f6480"}
N=230.0; ALPHA=12.34; rA=137.357; rB=126.75; Ry=133.0; EPD=25.9
a=np.radians(ALPHA)
bc=np.arcsin((0.75*1000)/(1e6/N)-np.sin(a))
RHO=(np.cos(a)+np.cos(bc))/(1.0/rA+1.0/rB)
des=R.Design("sim",N,ALPHA,rA,rB,Ry,RHO,EPD,BAND)
c=[0,0,0]; rBf=des.best_focus(c); C,ndet,edisp,ecross=des.detector(c,rBf)

# meridional trace (x=0 plane) -> (y,z) points, reusing the tracer's exact-sphere + local VLS
def meridional(lam,nrays=140):
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
    ic=nrays//2; chief=dout[ic]; beta=np.degrees(np.arctan2(chief[1],chief[2]))
    sel=ok&np.isfinite(H[:,0])
    yz=lambda M:np.stack([M[sel][:,2],M[sel][:,1]],axis=1)  # (z,y)
    Pyz=yz(P); Dyz=yz(dout); Dyz/=np.linalg.norm(Dyz,axis=1,keepdims=True)
    return np.array([des.A[2],des.A[1]]), Pyz, yz(H), beta, Dyz   # + diffracted dirs (z,y)

def focus_point(Pyz,Dyz):
    """least-squares point of closest approach of a fan of 2D rays (o_j, u_j) -> the wavelength focus."""
    M=np.zeros((2,2)); b=np.zeros(2)
    for o,u in zip(Pyz,Dyz):
        Pn=np.eye(2)-np.outer(u,u); M+=Pn; b+=Pn@o
    return np.linalg.solve(M,b)

# gather — true per-wavelength focus (line-of-closest-approach), not centroids on a fixed plane
Aslit,_,_,_,_=meridional(0.7)
foci=[]; betas={}
for l in BAND:
    _,P,H,beta,D=meridional(l); foci.append(focus_point(P,D)); betas[l]=beta
foci=np.array(foci)                          # (6,2) in (z,y) — the real focal curve
Ocen=np.array([0.,0.])                        # grating pole

# detector geometry from the ACTUAL focal locus
p=np.polyfit(foci[:,0],foci[:,1],1); zc=foci[:,0].mean(); yc=np.polyval(p,zc)
tdir=np.array([1,p[0]])/np.hypot(1,p[0]); half=29.18/2
e1=np.array([zc,yc])-half*tdir; e2=np.array([zc,yc])+half*tdir
spectrum_len=np.hypot(*(foci[-1]-foci[0]))    # 500->1000 nm run along the sensor
# tilt of detector vs the plane perpendicular to the central diffracted chief ray
cmid=np.array([zc,yc])/np.hypot(zc,yc); perp=np.array([-cmid[1],cmid[0]])
tilt=abs(np.degrees(np.arccos(np.clip(abs(tdir@perp),0,1))))
# detector angle from vertical (intuitive read)
det_from_vert=abs(np.degrees(np.arctan2(tdir[0],tdir[1])))

# ---------------- draw ----------------
fig,ax=plt.subplots(figsize=(15,9.2)); ax.set_facecolor("#0e1116")
CY="#39d98a"; DIM="#8fd3ff"; TXT="#dfe7ee"; GRY="#7f8c99"

# faint real ray fan (white in, coloured out)
_,P0,_,_,_=meridional(0.7,nrays=140)
for j in range(len(P0)): ax.plot([Aslit[0],P0[j,0]],[Aslit[1],P0[j,1]],color="#eaf2f8",lw=0.3,alpha=0.04,zorder=1)
for l in BAND:
    _,P,H,_,_=meridional(l,nrays=140); col=COL[int(l*1000)]
    for j in range(len(P)): ax.plot([P[j,0],H[j,0]],[P[j,1],H[j,1]],color=col,lw=0.3,alpha=0.06,zorder=1)

# grating surface (exact circle radius Ry) to scale
yy=np.linspace(-EPD/2-3,EPD/2+3,220); zz=Ry-np.sqrt(Ry*Ry-yy*yy)
ax.plot(zz,yy,color=CY,lw=3.2,zorder=6,solid_capstyle="round")
for yt in np.linspace(-EPD/2,EPD/2,24):
    zt=Ry-np.sqrt(Ry*Ry-yt*yt); ax.plot([zt,zt+1.0],[yt,yt],color=CY,lw=0.6,alpha=0.5,zorder=5)

# grating normal
ax.plot([0,zc*1.06],[0,0],ls=(0,(6,5)),color=GRY,lw=1.0,zorder=2)
ax.text(zc*1.06,0.6,"grating normal",color=GRY,fontsize=8.5,va="bottom",ha="right")

# entrance arm L_A (slit -> pole) + alpha
ax.plot([Aslit[0],0],[Aslit[1],0],color="#c9d2da",lw=1.3,zorder=4)
ax.add_patch(Arc((0,0),58,58,angle=0,theta1=0,theta2=ALPHA,color=DIM,lw=1.2))
ax.text(31,3.0,rf"$\alpha$ = {ALPHA:.2f}°",color=DIM,fontsize=10.5)
mid=0.5*Aslit; ax.text(mid[0],mid[1]+4,f"L$_A$ = {rA:.2f} mm",color="#c9d2da",fontsize=10,
        rotation=np.degrees(np.arctan2(Aslit[1],Aslit[0])),ha="center",va="bottom")

# slit: true 1 mm, perpendicular to the entrance arm
armdir=(-Aslit)/np.hypot(*Aslit); sperp=np.array([-armdir[1],armdir[0]])   # slit long axis
s1=Aslit+0.5*sperp; s2=Aslit-0.5*sperp                                     # 1 mm tall (true scale)
ax.plot([s1[0],s2[0]],[s1[1],s2[1]],color="#ffffff",lw=5,solid_capstyle="butt",zorder=8)
ax.annotate(f"Entrance slit — 20 µm × 1 mm\n⊥ entrance arm (plane tilted {ALPHA:.1f}° from vertical)\nfed by SMA905 + relay",
            xy=(Aslit[0],Aslit[1]),xytext=(Aslit[0]+18,Aslit[1]-10),color=TXT,fontsize=8.6,ha="right",va="top",
            arrowprops=dict(arrowstyle="->",color=GRY,lw=0.9))

# central diffracted arm L_H (pole -> detector centre)
ax.plot([0,zc],[0,yc],color="#c9d2da",lw=1.1,zorder=4)
mh=0.5*np.array([zc,yc]); ax.text(mh[0],mh[1]-4.5,f"L$_H$ = {rB:.2f} mm",color="#c9d2da",fontsize=9.5,
        rotation=np.degrees(np.arctan2(yc,zc)),ha="center",va="top")

# beta span arc between 500 and 1000 nm diffracted chief rays
r_arc=0.62*np.hypot(zc,yc)
b5=np.radians(betas[0.5]); b10=np.radians(betas[1.0])
ax.add_patch(Arc((0,0),2*r_arc,2*r_arc,angle=0,theta1=min(betas[0.5],betas[1.0]),theta2=max(betas[0.5],betas[1.0]),color="#c98f8f",lw=1.0))
ax.text(r_arc*np.cos((b5+b10)/2)+2,r_arc*np.sin((b5+b10)/2),
        rf"$\beta$: {betas[0.5]:.1f}° → {betas[1.0]:+.1f}°",color="#d7b2b2",fontsize=9,va="center")

# detector to scale on the real focal locus, at its true tilt
ax.plot([e1[0],e2[0]],[e1[1],e2[1]],color="#c9d2da",lw=6,solid_capstyle="butt",zorder=7)
for l in BAND:
    f=foci[BAND.index(l)]; ax.plot(f[0],f[1],'o',color=COL[int(l*1000)],ms=6,mec="k",mew=.4,zorder=9)
ax.text(foci[0][0]+2.0,foci[0][1]-1.5,"500 nm",color=COL[500],fontsize=8.5,ha="left",va="top")
ax.text(foci[-1][0]+2.0,foci[-1][1]+0.5,"1000 nm",color=COL[1000],fontsize=8.5,ha="left",va="bottom")
ax.annotate(f"Detector — flat, tilted {det_from_vert:.1f}° from vertical "
            f"({tilt:.1f}° off the central-ray normal)\n500–1000 nm spans {spectrum_len:.1f} mm on the array",
            xy=(e1[0],e1[1]),xytext=(zc-5,e1[1]-6),color=TXT,fontsize=8.6,ha="center",va="top",
            arrowprops=dict(arrowstyle="->",color=GRY,lw=0.9))

ax.set_xlim(-20,190); ax.set_ylim(-40,44)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("JASPER concave-VLS spectrometer core — TO-SCALE layout (true mm, equal aspect)   ·   slit & detector tilt measured off the ray trace",
             color="#eef3f7",fontsize=12.2,pad=14)
ax.text(0.5,-0.085,f"Concave VLS grating 230 l/mm, R = {Ry:.0f} mm, m = 1.  Finite-arm Shimadzu-style self-imaging mount: "
        f"L$_A$ = {rA:.1f} mm above the normal, L$_H$ = {rB:.1f} mm below, plate scale 34.2 nm/mm.  "
        f"Meridional plane — sphere and toroid coincide here; the toroid differs only in astigmatic line height.",
        transform=ax.transAxes,color="#9fb0bf",fontsize=8.5,ha="center",va="top",wrap=True)
fig.tight_layout(); fig.savefig("jasper_layout_toscale.png",dpi=155,bbox_inches="tight",facecolor="#0e1116")
print(f"wrote jasper_layout_toscale.png | LA={rA:.2f} LH={rB:.2f} | beta {betas[0.5]:.2f}->{betas[1.0]:.2f} | "
      f"det tilt {det_from_vert:.1f} from vert, {tilt:.1f} off central-normal | spectrum {spectrum_len:.2f} mm")
