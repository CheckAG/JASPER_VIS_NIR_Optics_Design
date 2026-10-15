import numpy as np, json, vls_raytrace as R, validate_p0400 as V
from scipy.optimize import minimize
N=400.; A=20.; rA=96.9; RS=V.RSUB; EPD=18.3; BAND=V.BAND; CS=V.CS
a=np.radians(A); bref=np.radians(V.PBETA); D=V.PDIST
beta=lambda l,nu=N: np.arcsin(l*1e-3*nu-np.sin(a))
rp=lambda l,nu=N: D/np.cos(beta(l,nu)-bref)

# ---------- Q1: where is the 400? test the alternatives against the datasheet ----------
cfit=[0.2528,-0.01713,-3.0e-7]
nu_at=lambda y: N+cfit[0]*y+cfit[1]*y*y+cfit[2]*y**3
print("Q1  hypothesis test (datasheet: beta -10.49/+5.62 deg, L 27.7 mm, 25.3 nm/mm)")
for tag,nupole in (("400 at the POLE",400.0),("400 at edge y=-20",400+(400-nu_at(-20))),("400 at edge y=+20",400+(400-nu_at(20)))):
    b1,b2=np.degrees(beta(0.4,nupole)),np.degrees(beta(1.1,nupole))
    p1=rp(0.4,nupole)*np.array([np.sin(beta(0.4,nupole)),np.cos(beta(0.4,nupole))]); p2=rp(1.1,nupole)*np.array([np.sin(beta(1.1,nupole)),np.cos(beta(1.1,nupole))])
    L=np.linalg.norm(p2-p1); disp=700/L
    print(f"  {tag:<20} pole density {nupole:6.2f} -> beta {b1:+.2f}/{b2:+.2f} deg, L {L:.2f} mm, {disp:.2f} nm/mm")
dbdnu=np.degrees(1.1e-3/np.cos(beta(1.1)))
print(f"  sensitivity: 1 l/mm at the pole moves beta(1100 nm) by {dbdnu:.3f} deg -> datasheet angles (0.01 deg) pin the pole to +/-{0.01/dbdnu:.2f} l/mm")
# orientation: blank mounted 180 deg rotated (y -> -y): odd terms flip
def dl_mean(c):
    des=R.Design("p",N,A,rA,rA,RS,RS,EPD,BAND); rows,_=V.resolution(des,c,EPD); return np.array([r['dl'] for r in rows])
ok=dl_mean(cfit); fl=dl_mean([-cfit[0],cfit[1],-cfit[2]])
print(f"  orientation: correct {ok.mean():.2f} nm mean | blank rotated 180 deg {fl.mean():.2f} nm mean (worst {fl.max():.1f})")

# ---------- Q2: decompose 0.62 (LPF idealised) vs 2.02 (ray trace, one chirp) ----------
print("\nQ2  per-wavelength resolution (nm), f/5.3, datasheet plane, R=%.2f"%RS)
des=R.Design("p",N,A,rA,rA,RS,RS,EPD,BAND)
def dl_one(l,c):
    Cs=CS; U,_=des.trace(l,c,npup=2500,Csys=Cs); rms=(U-U.mean()).std()*1e3
    b=beta(l); r=rp(l); plate=(1e6/N)*np.cos(b)*np.cos(b-bref)/r   # nm/mm along the tilted plane, Shimadzu eq (24)
    diff=1.028*l*r/EPD; return np.hypot(np.hypot(20,2.355*rms),np.hypot(8,diff))/1e3*plate, rms
# (a) tracer, idealised: best chirp for EACH wavelength separately
sc=np.array([0.1,0.01,1e-4]); per=[]
for l in BAND:
    f=lambda z:(lambda U:(U-U.mean()).std())(des.trace(l,list(z*sc),npup=400,Csys=CS)[0])
    z=min((minimize(f,s,method="Nelder-Mead",options=dict(xatol=1e-5,fatol=1e-7,maxiter=1500)) for s in ([2.5,-1.7,0],[1,-1,0])),key=lambda r:r.fun).x
    per.append(dl_one(l,list(z*sc))[0])
# (b) tracer, real: one chirp for the whole band
one=[dl_one(l,cfit)[0] for l in BAND]
# (c) analytic LPF, idealised (F20 = 0 at every lambda, coma only) and (d) analytic with the SAME single c1
def F20(b,r): return np.cos(a)**2/rA-np.cos(a)/RS+np.cos(b)**2/r-np.cos(b)/RS
def F30(b,r): return (np.sin(a)/rA)*(np.cos(a)**2/rA-np.cos(a)/RS)+(np.sin(b)/r)*(np.cos(b)**2/r-np.cos(b)/RS)
w=np.linspace(-EPD/2,EPD/2,201); lpf_id=[]; lpf_one=[]
c1n=[F20(beta(l),rp(l))/(l*1e-3) for l in BAND]
for l,cn in zip(BAND,c1n):
    b=beta(l); r=rp(l); plate=(1e6/N)*np.cos(b)*np.cos(b-bref)/r; diff=1.028*l*r/EPD
    def dlam(F2):
        dy=(r/np.cos(b))*(F2*w+1.5*F30(b,r)*w**2); dy-=dy.mean(); rms=np.sqrt(np.mean(dy**2))*1e3
        return np.hypot(np.hypot(20,2.355*rms),np.hypot(8,diff))/1e3*plate
    lpf_id.append(dlam(0.0)); lpf_one.append(dlam(l*1e-3*(cn+abs(cfit[0]))))   # analytic sign convention: c1 = -|c1_trace|
lam=[int(l*1e3) for l in BAND]
for name,v in (("LPF idealised (F20=0 every lambda)",lpf_id),("TRACE idealised (best chirp per lambda)",per),
               ("LPF with the one fitted c1",lpf_one),("TRACE with the one fitted chirp",one)):
    print(f"  {name:<40}"+" ".join(f"{x:5.2f}" for x in v)+f"   mean {np.mean(v):.2f}")
json.dump(dict(lam=lam,lpf_id=lpf_id,trace_id=per,lpf_one=lpf_one,trace_one=one,ok=list(ok),flipped=list(fl)),open("p0400_verify.json","w"),indent=1,default=float)
