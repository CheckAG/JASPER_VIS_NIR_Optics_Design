"""VLS / holographic grating on a CURVED substrate for Optiland.
Implements the capability Optiland lacks: variable-line-space grooves on a spherical/conic
substrate. Done by swapping a grating surface's interaction model for one that computes a
POSITION-DEPENDENT local groove spacing d(y)=1/nu(y) along the dispersion axis (y).

  nu(y) = nu0 + c1*y + c2*y^2 + c3*y^3     grooves per MICROMETRE   (y in mm, local surface coord)

UNITS NOTE: Optiland's grating_period is in micrometres while positions/radii are in mm, so
nu0 is in grooves/um:  600 l/mm = 0.6 /um,  300 l/mm = 0.3 /um.  Constant c's -> plain grating.
"""
import numpy as np, optiland.backend as be
from optiland.interactions.diffractive_model import DiffractiveInteractionModel

class VLSDiffractiveModel(DiffractiveInteractionModel):
    interaction_type = "vls_diffractive"
    def set_vls(self, nu0, c1=0.0, c2=0.0, c3=0.0):
        self.nu0=nu0; self.c1=c1; self.c2=c2; self.c3=c3   # nu0 in grooves/um
    def interact_real_rays(self, rays):
        nx,ny,nz = self.geometry.surface_normal(rays)
        n1 = self.material_pre.n(rays.w); n2 = self.material_post.n(rays.w)
        fx,fy,fz = self.geometry.grating_vector(rays)
        y = rays.y
        nu = self.nu0 + self.c1*y + self.c2*y*y + self.c3*y*y*y   # local grooves/um along dispersion axis
        period_local = 1.0/nu                                    # um (matches Optiland grating_period units)
        pp = period_local / be.sqrt(fx**2 + fy**2)
        m = self.geometry.grating_order
        rays.gratingdiffract(nx,ny,nz,fx,fy,fz,m,pp,n1,n2,self.is_reflective)
        rays = self._apply_coating_and_bsdf(rays, nx,ny,nz)
        return rays

def make_vls(surface, nu0, c1=0.0,c2=0.0,c3=0.0):
    """Convert an existing grating surface's interaction model into a VLS one in place."""
    im=surface.interaction_model
    im.__class__=VLSDiffractiveModel
    im.set_vls(nu0,c1,c2,c3)
    return surface

if __name__=="__main__":
    # VALIDATION: constant nu (c=0) must reproduce the plain grating equation exactly.
    from optiland.optic import Optic
    nu0=0.600    # grooves per MICROMETRE == 600 l/mm
    o=Optic()
    o.surfaces.add(index=0,radius=be.inf,thickness=be.inf)
    o.surfaces.add(index=1,surface_type="grating",material="mirror",radius=be.inf,thickness=-50.0,
                   is_stop=True,grating_order=1,grating_period=1.0/nu0,groove_orientation_angle=0.0)
    o.surfaces.add(index=2)
    o.set_aperture(aperture_type="EPD",value=4.0)
    o.fields.set_type(field_type="angle"); o.fields.add(x=0,y=0)
    for w in (0.5,0.7,1.0): o.wavelengths.add(value=w)
    make_vls(o.surfaces.surfaces[1], nu0=nu0)   # constant -> should match grating eq
    g=lambda a: float(np.asarray(be.to_numpy(a)).ravel()[0])
    print("VLS constant-nu check vs grating equation (600 l/mm, normal incidence):")
    mx=0.0
    for w in (0.5,0.7,1.0):
        r=o.trace_generic(Hx=0,Hy=0,Px=0.0,Py=0.0,wavelength=w)
        M=abs(g(r.M)); pred=w*nu0        # lam(um) * nu0(1/um) = lam/d
        mx=max(mx,abs(M-pred))
        print(f"  {w:.2f}um  |M|={M:.5f}  lam/d={pred:.5f}")
    print(f"  max err={mx:.2e}  ({'PASS' if mx<1e-9 else 'CHECK'})")
