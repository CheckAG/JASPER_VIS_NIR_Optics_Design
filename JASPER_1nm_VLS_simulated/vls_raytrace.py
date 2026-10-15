#!/usr/bin/env python3
"""Standalone geometric ray tracer for a concave VLS grating in a finite-arm mount — NO Optiland.

Traces real rays: slit point-source -> intersect the toroidal grating surface -> apply the LOCAL
variable-line-spacing vector grating equation at the hit point -> propagate to a flat detector
placed on the diffracted chief ray -> collect (dispersion, cross-dispersion) hits. From the hits it
computes per-wavelength spot size, spectral resolution, plate scale and astigmatic line height.

Purpose: give Optiland-independent, geometry-based spot diagrams WITH variable line spacing, and
cross-validate the analytic light-path-function (LPF) numbers. Two modes per design:
  - CLASSICAL grooves (constant density): should reproduce the LPF 'classical' (uncorrected) numbers
    with NO assumptions -> a clean check of the ray tracer against the analytic theory.
  - FLAT-FIELD: fit the VLS coefficients (c1,c2,c3) to minimise dispersion RMS -> best-case.

Coordinates (mm): grating pole at origin O; normal +z; dispersion plane = y-z; grooves along x.
Slit at angle alpha from the normal, at +y. All lengths mm, wavelengths mm, groove density /mm.
"""
import numpy as np
from scipy.optimize import minimize

# ---------------- grating surface (paraxial toroid: matches the LPF surface) ----------------
def sag(x, y, Ry, rho):            # z of the surface at (x,y); concave toward +z
    return x*x/(2*rho) + y*y/(2*Ry)
def surf_normal(x, y, Ry, rho):    # unit normal from gradient of (z - sag)
    n = np.array([-x/rho, -y/Ry, 1.0]); return n/np.linalg.norm(n)
def intersect(A, d, Ry, rho):      # ray A + t d meets the surface; Newton from the z=0 crossing
    t = -A[2]/d[2]
    for _ in range(60):
        P = A + t*d
        S = P[2] - sag(P[0], P[1], Ry, rho)
        dS = d[2] - (P[0]/rho*d[0] + P[1]/Ry*d[1])
        step = S/dS; t -= step
        if abs(step) < 1e-12: break
    return A + t*d

# ---------------- local VLS vector grating equation (reflection) ----------------
def diffract(d_in, P, N, nu_local, m, lam):
    # grating vector g: in-surface unit vector perpendicular to the grooves (grooves along x)
    yhat = np.array([0.0,1.0,0.0])
    g = yhat - np.dot(yhat, N)*N; g = g/np.linalg.norm(g)
    din_tan = d_in - np.dot(d_in, N)*N
    dout_tan = din_tan + m*lam*nu_local*g            # add grating momentum along g
    n2 = 1.0 - np.dot(dout_tan, dout_tan)
    if n2 < 0: return None                            # evanescent -> lost
    # reflection: normal component opposite sign to the incident normal component
    sign = -np.sign(np.dot(d_in, N))
    return dout_tan + sign*np.sqrt(n2)*N

def nu_of(y, nu0, c1, c2, c3):     # local groove density (/mm) vs dispersion coordinate y (mm)
    return nu0 + c1*y + c2*y*y + c3*y*y*y

# ---------------- one design ----------------
class Design:
    def __init__(s, name, N, alpha_deg, rA, rB, Ry, rho, EPD, band, m=1):
        s.name=name; s.nu0=N; s.a=np.radians(alpha_deg); s.rA=rA; s.rB=rB
        s.Ry=Ry; s.rho=rho; s.EPD=EPD; s.band=np.array(band); s.m=m
        s.A=np.array([0.0, rA*np.sin(s.a), rA*np.cos(s.a)])         # slit point (mm)
    def chief_out(s, lam, c):
        d_in=(np.zeros(3)-s.A); d_in/=np.linalg.norm(d_in)         # slit -> pole
        N=np.array([0.,0.,1.])
        return diffract(d_in, np.zeros(3), N, s.nu0, s.m, lam*1e-3)  # lam um->mm; pole nu=nu0
    def detector(s, c, rB=None):
        rB=s.rB if rB is None else rB
        lam_c=s.band[len(s.band)//2]
        dc=s.chief_out(lam_c, c)
        C=rB*dc                                                    # detector centre on chief ray at rB
        ndet=-dc                                                   # faces the beam
        ecross=np.array([1.,0.,0.])                               # grooves axis (x) in detector plane
        edisp=np.cross(ndet, ecross); edisp/=np.linalg.norm(edisp)
        return C, ndet, edisp, ecross
    def best_focus(s, c):                                          # detector distance that best-focuses the band
        def rms_at(rB):
            Csys=s.detector(c, rB); out=[]
            for lam in s.band:
                U,_=s.trace(lam,c,npup=250,Csys=Csys)
                if len(U)<40: return 1e6
                out.append((U-U.mean()).std())
            return float(np.mean(out))
        from scipy.optimize import minimize_scalar
        r=minimize_scalar(rms_at, bounds=(0.9*s.rB,1.1*s.rB), method="bounded")
        return r.x
    def trace(s, lam, c, npup=900, Csys=None):
        if Csys is None: Csys=s.detector(c)
        C,ndet,edisp,ecross=Csys
        Ry,rho,A,lam_mm=s.Ry,s.rho,s.A,lam*1e-3
        exact_sphere = abs(rho-Ry)<1e-9
        # sample pupil = grating aperture disk (radius EPD/2) in x,y  (deterministic per npup)
        rng=np.random.default_rng(7)
        th=rng.uniform(0,2*np.pi,npup); rr=np.sqrt(rng.uniform(0,1,npup))*s.EPD/2
        xs=rr*np.cos(th); ys=rr*np.sin(th)
        # aim points on the surface (for the ray direction); sag is a fine aim target
        Pt=np.stack([xs,ys,xs*xs/(2*rho)+ys*ys/(2*Ry)],axis=1)     # (N,3)
        d=Pt-A; d/=np.linalg.norm(d,axis=1,keepdims=True)
        if exact_sphere:
            R=Ry; Oc=np.array([0.,0.,R])                           # sphere centre (concave toward +z)
            e=A-Oc; ed=d@e; disc=ed*ed-(e@e-R*R)
            t=-ed+np.sqrt(np.clip(disc,0,None))                    # slit is inside the sphere -> far root
            P=A+t[:,None]*d
            N=(Oc-P)/R                                             # unit normal, +z at pole
        else:
            # true toroid: meridional circle radius Ry (y-z), sagittal radius rho (x); Newton on implicit F
            t=-A[2]/d[:,2]
            for _ in range(60):
                P=A+t[:,None]*d
                # z of a toroid: sphere-exact in y, circular in x about the meridian
                zy=Ry-np.sqrt(np.clip(Ry*Ry-P[:,1]**2,1e-9,None))   # meridional sag (exact circle)
                F=(rho-(P[:,2]-zy))**2 + P[:,0]**2 - rho*rho        # sagittal circle of radius rho
                dzy=P[:,1]/np.sqrt(np.clip(Ry*Ry-P[:,1]**2,1e-9,None))
                dF=2*(rho-(P[:,2]-zy))*(-(d[:,2]-dzy*d[:,1])) + 2*P[:,0]*d[:,0]
                step=F/dF; t=t-step
            P=A+t[:,None]*d
            zy=Ry-np.sqrt(np.clip(Ry*Ry-P[:,1]**2,1e-9,None))
            dzy=P[:,1]/np.sqrt(np.clip(Ry*Ry-P[:,1]**2,1e-9,None))
            # normal from grad of F(x,y,z): dF/dx=2x, dF/dz=-2(rho-(z-zy)), dF/dy=2(rho-(z-zy))*dzy
            k=(rho-(P[:,2]-zy))
            N=np.stack([-P[:,0], -k*dzy, k],axis=1)                # oriented +z
            N/=np.linalg.norm(N,axis=1,keepdims=True)
        # local VLS grating vector g (perp to grooves, in surface) and density
        yhat=np.array([0.,1.,0.])
        g=yhat-(N[:,1:2])*N; g/=np.linalg.norm(g,axis=1,keepdims=True)
        nu=nu_of(P[:,1],s.nu0,*c)
        din_tan=d-(np.sum(d*N,axis=1,keepdims=True))*N
        dout_tan=din_tan+(s.m*lam_mm*nu)[:,None]*g
        n2=1.0-np.sum(dout_tan*dout_tan,axis=1)
        ok=n2>=0
        sign=-np.sign(np.sum(d*N,axis=1))
        dout=dout_tan+(sign*np.sqrt(np.clip(n2,0,None)))[:,None]*N
        # intersect detector plane {X: (X-C).ndet=0}
        denom=dout@ndet
        tt=((C-P)@ndet)/denom
        H=P+tt[:,None]*dout
        U=(H-C)@edisp; V=(H-C)@ecross
        m=ok&np.isfinite(U)&np.isfinite(V)
        return U[m],V[m]

SLIT=20.0; PIX=8.0  # um

def evaluate(des, c, label, rB=None):
    if rB is None: rB=des.best_focus(c)
    Csys=des.detector(c, rB); centroids=[]; rows=[]
    for lam in des.band:
        U,V=des.trace(lam,c,Csys=Csys)
        cu=U.mean(); centroids.append(cu)
    # plate scale from centroid vs lambda (nm per mm); local slope
    lam_nm=des.band*1e3                                            # um -> nm
    for i,lam in enumerate(des.band):
        U,V=des.trace(lam,c,Csys=Csys)
        rms_disp=(U-U.mean()).std()*1000.0                        # mm -> um
        line=(V.max()-V.min())*1.0                                 # mm (cross-dispersion span)
        a=max(0,i-1); b=min(len(des.band)-1,i+1)
        plate=abs((lam_nm[b]-lam_nm[a])/((centroids[b]-centroids[a])))  # nm/mm
        diff=1.028*lam*(des.rB/des.EPD)                            # diffraction FWHM (um), lam in um, f/# = rB/EPD
        dlam=np.hypot(np.hypot(SLIT,2.355*rms_disp),np.hypot(PIX,diff))/1000.0*plate
        rows.append(dict(lam=int(lam*1e3),rms=rms_disp,line=line,plate=plate,dlam=dlam))
    span=abs(centroids[0]-centroids[-1])
    mean=np.mean([r['dlam'] for r in rows]); worst=max(r['dlam'] for r in rows)
    aline=np.mean([r['line'] for r in rows])
    print(f"\n[{des.name}] {label}")
    print(f"  detector focus r'={rB:.2f} mm | plate {rows[0]['plate']:.1f}-{rows[-1]['plate']:.1f} nm/mm | span {span:.1f} mm | "
          f"Dlam mean {mean:.2f} worst {worst:.2f} nm | astig line {aline:.2f} mm")
    for r in rows: print(f"    {r['lam']:4d} nm  RMSdisp {r['rms']:5.1f} um  line {r['line']:4.2f} mm  dlam {r['dlam']:.2f} nm")
    return dict(rows=rows,mean=mean,worst=worst,span=span,aline=aline)

def fit_vls(des):
    def obj(c):
        rB=des.best_focus(c); Csys=des.detector(c, rB); out=[]
        for lam in des.band:
            U,_=des.trace(lam,c,npup=250,Csys=Csys)
            if len(U)<40: return 1e6
            out.append((U-U.mean()).std()*1000)
        out=np.array(out); return out.mean()+0.6*out.max()
    best=None
    for seed in ([0,0,0],[1e-3,1e-5,0]):
        r=minimize(obj,seed,method="Nelder-Mead",options=dict(xatol=1e-6,fatol=0.1,maxiter=250))
        if best is None or r.fun<best.fun: best=r
    return best.x

if __name__=="__main__":
    band=[0.5,0.6,0.7,0.8,0.9,1.0]  # um
    # Shimadzu P0400-01 (real, spherical): 400 l/mm, alpha=20, arms 96.9, R~100.3
    p04=Design("Shimadzu P0400-01",400.0,20.0,96.9,96.9,100.29,100.29,18.3,band)
    evaluate(p04,[0,0,0],"CLASSICAL grooves  [detector here is perpendicular to the mid-band ray, NOT the tilted datasheet plane -- not comparable to the LPF number; see validate_models.py]")
    c=fit_vls(p04); print(f"  fitted VLS c={c}")
    evaluate(p04,c,"FLAT-FIELD (fitted VLS)  [LPF ~0.60 nm is an idealisation (defocus nulled at every lambda); see validate_models.py]")
