import json, math, os
HERE=os.path.dirname(os.path.abspath(__file__)); DESIGN=HERE   # release: one directory
S=json.load(open(os.path.join(HERE,"spots.json")))  # produced by analyze_spots.py
WAVES=[500,600,700,800,900,1000]   # 750 nm not in the MED/RAY for this build
COL={500:"#1aa5cf",600:"#e8912a",700:"#e34b3a",800:"#b23a44",900:"#8f4f63",1000:"#6f6480"}
NA=4.875      # image-space f/# = 78/16 (f/4.9)
SLIT_UM=30.0  # entrance slit
PIX_UM=8.0    # TCD1304: 3648 px x 8 um = 29.1 mm
ALPHA=35.64   # grating incidence, deg
GROOVES=600.0
DET_MM=29.1   # TCD1304 active length

def plate_scale(w):
    """Local nm/mm from the traced centroids -- nonlinear, so never use a band average."""
    i=WAVES.index(w); a,b=WAVES[max(0,i-1)],WAVES[min(len(WAVES)-1,i+1)]
    return abs((b-a)/(S[f"{b}_ctr"]["cx"]-S[f"{a}_ctr"]["cx"]))

def slit_image_um(w):
    """Slit width at the detector: unit magnification, narrowed by the anamorphic factor."""
    beta=math.degrees(math.asin(w*1e-6*GROOVES-math.sin(math.radians(ALPHA))))
    return SLIT_UM*math.cos(math.radians(ALPHA))/math.cos(math.radians(beta))
def esc(s): return str(s)
# gather center-field data
rows=[]
for w in WAVES:
    d=S[f"{w}_ctr"]
    airy=1.22*(w/1000.0)*NA           # um radius
    # delivered resolution: slit image (+) geometric blur (+) pixel, in quadrature,
    # converted with the LOCAL plate scale (CLAUDE.md "Physics conventions")
    tot=math.hypot(math.hypot(slit_image_um(w),2.355*d["rmsx"]),PIX_UM)
    dlam=tot/1000.0*plate_scale(w)
    rows.append(dict(w=w,n=d["n"],cx=d["cx"],rmsr=d["rmsr"],rmsx=d["rmsx"],rmsy=d["rmsy"],
                     geo=d["geo"],airy=airy,dlam=dlam,pts=d["pts"],col=COL[w]))
maxgeo=max(r["geo"] for r in rows)      # 389.7
HALF=430.0                              # common panel half-scale (um)

# ---------- detector strip ----------
W=1040; SL=60; SR=60; y0=46; sh=30
plotw=W-SL-SR
def mm2x(mm): return SL+(mm+DET_MM/2)/DET_MM*plotw
strip=[f'<svg viewBox="0 0 {W} 118" width="100%" role="img" aria-label="Spectrum on the detector">']
strip.append(f'<rect x="{SL}" y="{y0}" width="{plotw}" height="{sh}" rx="4" class="sensor"/>')
# pixel ticks every ~2mm
for mm in range(-14,15,2):
    x=mm2x(mm); strip.append(f'<line x1="{x:.1f}" y1="{y0+sh}" x2="{x:.1f}" y2="{y0+sh+4}" class="axtick"/>')
    strip.append(f'<text x="{x:.1f}" y="{y0+sh+16}" class="axlab" text-anchor="middle">{mm}</text>')
strip.append(f'<text x="{SL+plotw/2:.1f}" y="{y0+sh+30}" class="axtitle" text-anchor="middle">detector position (mm) — 29.1 mm / 3648 px</text>')
for r in rows:
    x=mm2x(r["cx"])
    strip.append(f'<line x1="{x:.1f}" y1="{y0-8}" x2="{x:.1f}" y2="{y0+sh}" stroke="{r["col"]}" stroke-width="2.4"/>')
    strip.append(f'<circle cx="{x:.1f}" cy="{y0-8}" r="4.5" fill="{r["col"]}"/>')
    strip.append(f'<text x="{x:.1f}" y="{y0-16}" class="wlab" text-anchor="middle" fill="{r["col"]}">{r["w"]}</text>')
strip.append('</svg>')
strip="".join(strip)

# ---------- RMS bar chart ----------
BW=1040; bl=64; br=24; bt=24; bb=54; ph=250
pw=BW-bl-br; pph=ph-bt-bb
ymax=60
def by(v): return bt+pph-(v/ymax)*pph
bar=[f'<svg viewBox="0 0 {BW} {ph}" width="100%" role="img" aria-label="RMSx spot radius by wavelength">']
for gv in [0,10,20,30,40,50]:
    yy=by(gv); bar.append(f'<line x1="{bl}" y1="{yy:.1f}" x2="{BW-br}" y2="{yy:.1f}" class="grid"/>')
    bar.append(f'<text x="{bl-10}" y="{yy+4:.1f}" class="axlab" text-anchor="end">{gv}</text>')
bar.append(f'<text x="16" y="{bt+pph/2:.1f}" class="axtitle" text-anchor="middle" transform="rotate(-90 16 {bt+pph/2:.1f})">RMSx — dispersion axis (µm)</text>')
# pixel reference line 14um
yp=by(8); bar.append(f'<line x1="{bl}" y1="{yp:.1f}" x2="{BW-br}" y2="{yp:.1f}" class="refline"/>')
bar.append(f'<text x="{BW-br-4}" y="{yp-5:.1f}" class="reflab" text-anchor="end">8 µm detector pixel (TCD1304)</text>')
n=len(rows); slot=pw/n; bwid=slot*0.52
for i,r in enumerate(rows):
    cx=bl+slot*(i+0.5); x=cx-bwid/2; h=(r["rmsx"]/ymax)*pph; yy=by(r["rmsx"])
    bar.append(f'<rect x="{x:.1f}" y="{yy:.1f}" width="{bwid:.1f}" height="{h:.1f}" rx="3" fill="{r["col"]}"><title>{r["w"]} nm — RMSx {r["rmsx"]:.0f} µm, RMSy {r["rmsy"]:.0f} µm, Δλ {r["dlam"]:.1f} nm, n={r["n"]}/169</title></rect>')
    bar.append(f'<text x="{cx:.1f}" y="{yy-6:.1f}" class="barval" text-anchor="middle" fill="{r["col"]}">{r["rmsx"]:.0f}</text>')
    bar.append(f'<text x="{cx:.1f}" y="{bt+pph+18:.1f}" class="axlab" text-anchor="middle">{r["w"]}</text>')
bar.append(f'<text x="{bl+pw/2:.1f}" y="{ph-6:.1f}" class="axtitle" text-anchor="middle">wavelength (nm)</text>')
bar.append('</svg>')
bar="".join(bar)

# ---------- spot panels ----------
PS=176; pad=20; plot=PS-2*pad; cxp=PS/2; sc=(plot/2)/HALF
def panel(r):
    s=[f'<svg viewBox="0 0 {PS} {PS}" width="100%" role="img" aria-label="{r["w"]} nm spot">']
    s.append(f'<rect x="{pad}" y="{pad}" width="{plot}" height="{plot}" class="pbox"/>')
    s.append(f'<line x1="{cxp}" y1="{pad}" x2="{cxp}" y2="{PS-pad}" class="pcross"/>')
    s.append(f'<line x1="{pad}" y1="{cxp}" x2="{PS-pad}" y2="{cxp}" class="pcross"/>')
    rmsR=r["rmsr"]*sc
    s.append(f'<circle cx="{cxp}" cy="{cxp}" r="{rmsR:.1f}" class="rmsring" stroke="{r["col"]}"/>')
    for px,py in r["pts"]:
        x=cxp+px*sc; y=cxp-py*sc
        s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.5" fill="{r["col"]}" fill-opacity="0.85"/>')
    ar=max(r["airy"]*sc,0.8)
    s.append(f'<circle cx="{cxp}" cy="{cxp}" r="{ar:.2f}" class="airy"/>')
    # scale bar 100um
    sb=100*sc; sx=PS-pad-sb-2; sy=PS-pad-6
    s.append(f'<line x1="{sx:.1f}" y1="{sy:.1f}" x2="{sx+sb:.1f}" y2="{sy:.1f}" class="scalebar"/>')
    s.append(f'<text x="{sx+sb/2:.1f}" y="{sy-4:.1f}" class="scalelab" text-anchor="middle">100 µm</text>')
    s.append('</svg>')
    return "".join(s)
panels=[]
for r in rows:
    # regime by delivered resolution, not by raw spot size
    tag = ('<span class="chip good">near best focus</span>' if r["dlam"]<2.5
           else '<span class="chip crit">residual defocus</span>' if r["dlam"]>4.0
           else '<span class="chip warn">part defocused</span>')
    panels.append(f'''<figure class="panel">
      <div class="ptitle"><span class="pw" style="color:{r["col"]}">{r["w"]} nm</span>{tag}</div>
      {panel(r)}
      <figcaption>
        <div><span>RMSx</span><b>{r["rmsx"]:.0f} µm</b></div>
        <div><span>GEO</span><b>{r["geo"]:.0f} µm</b></div>
        <div><span>rays</span><b>{r["n"]}/169</b></div>
        <div><span>Δλ</span><b>{r["dlam"]:.1f} nm</b></div>
      </figcaption>
    </figure>''')
panels="".join(panels)

best=min(rows,key=lambda r:r["dlam"])
span=rows[0]["cx"]-rows[-1]["cx"]
minn=min(r["n"] for r in rows); maxn=max(r["n"] for r in rows)
surv=sum(r["n"] for r in rows)
html=f'''<meta charset="utf-8"><style>
:root{{
  --bg:#f6f8f9; --surface:#ffffff; --panel:#fbfcfc; --ink:#12181a; --ink2:#4a565a;
  --muted:#8a959a; --line:#e3e8ea; --hair:#eef2f3; --accent:#0f9a8c; --sensor:#dfe6e8;
  --refline:#c08a2e;
  --font:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  --mono:ui-monospace,"SF Mono","Cascadia Mono",Menlo,Consolas,monospace;
}}
@media (prefers-color-scheme:dark){{:root{{
  --bg:#0e1315; --surface:#141a1c; --panel:#161d1f; --ink:#eef3f4; --ink2:#a7b3b6;
  --muted:#71807f; --line:#232c2e; --hair:#1c2426; --accent:#1fc0ae; --sensor:#232c2e;
  --refline:#d9a441;
}}}}
:root[data-theme="light"]{{
  --bg:#f6f8f9; --surface:#ffffff; --panel:#fbfcfc; --ink:#12181a; --ink2:#4a565a;
  --muted:#8a959a; --line:#e3e8ea; --hair:#eef2f3; --accent:#0f9a8c; --sensor:#dfe6e8; --refline:#c08a2e;
}}
:root[data-theme="dark"]{{
  --bg:#0e1315; --surface:#141a1c; --panel:#161d1f; --ink:#eef3f4; --ink2:#a7b3b6;
  --muted:#71807f; --line:#232c2e; --hair:#1c2426; --accent:#1fc0ae; --sensor:#232c2e; --refline:#d9a441;
}}
*{{box-sizing:border-box}}
.wrap{{background:var(--bg);color:var(--ink);font-family:var(--font);
  -webkit-font-smoothing:antialiased;line-height:1.5;min-height:100vh;padding:34px 22px 60px}}
.inner{{max-width:1100px;margin:0 auto}}
.top{{display:flex;justify-content:space-between;align-items:flex-start;gap:20px;flex-wrap:wrap}}
.eyebrow{{font-family:var(--mono);font-size:11.5px;letter-spacing:.16em;text-transform:uppercase;
  color:var(--accent);margin:0 0 8px}}
h1{{font-size:29px;line-height:1.12;margin:0;letter-spacing:-.01em;text-wrap:balance;max-width:20ch}}
.sub{{color:var(--ink2);margin:10px 0 0;max-width:64ch;font-size:14.5px}}
.tog{{font-family:var(--mono);font-size:11px;letter-spacing:.08em;text-transform:uppercase;
  background:var(--surface);color:var(--ink2);border:1px solid var(--line);border-radius:7px;
  padding:8px 12px;cursor:pointer}}
.tog:hover{{color:var(--ink);border-color:var(--accent)}}
.tiles{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:26px 0 10px}}
.tile{{background:var(--surface);border:1px solid var(--line);border-radius:11px;padding:15px 16px}}
.tile .k{{font-family:var(--mono);font-size:10.5px;letter-spacing:.1em;text-transform:uppercase;color:var(--muted)}}
.tile .v{{font-size:25px;font-weight:600;margin-top:5px;letter-spacing:-.01em;font-variant-numeric:tabular-nums}}
.tile .v small{{font-size:14px;font-weight:500;color:var(--ink2)}}
.card{{background:var(--surface);border:1px solid var(--line);border-radius:13px;padding:20px 22px;margin-top:20px}}
.card h2{{font-size:16px;margin:0 0 3px;letter-spacing:-.005em}}
.card .cap{{color:var(--ink2);font-size:13px;margin:0 0 16px;max-width:70ch}}
.sensor{{fill:var(--sensor);stroke:var(--line)}}
.axtick{{stroke:var(--muted)}} .grid{{stroke:var(--hair)}}
.axlab{{font-family:var(--mono);font-size:10.5px;fill:var(--muted)}}
.axtitle{{font-family:var(--mono);font-size:10.5px;fill:var(--ink2);letter-spacing:.04em}}
.wlab{{font-family:var(--mono);font-size:11px;font-weight:600}}
.barval{{font-family:var(--mono);font-size:11px;font-weight:600;font-variant-numeric:tabular-nums}}
.refline{{stroke:var(--refline);stroke-width:1.3;stroke-dasharray:5 4}}
.reflab{{font-family:var(--mono);font-size:10px;fill:var(--refline)}}
.matrix{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px}}
.panel{{margin:0;background:var(--panel);border:1px solid var(--line);border-radius:11px;padding:12px 12px 10px}}
.ptitle{{display:flex;justify-content:space-between;align-items:flex-start;margin-bottom:4px;min-height:30px}}
.pw{{font-family:var(--mono);font-weight:600;font-size:13px;white-space:nowrap}}
.chip{{font-family:var(--mono);font-size:9px;letter-spacing:.04em;text-transform:uppercase;line-height:1.25;
  width:62px;text-align:center;white-space:normal;
  padding:2px 4px;border-radius:5px;border:1px solid transparent}}
.chip.good{{color:#0ca30c;border-color:#0ca30c55}} .chip.warn{{color:#c98500;border-color:#c9850055}}
.chip.crit{{color:#d64545;border-color:#d6454555}}
.pbox{{fill:none;stroke:var(--line)}}
.pcross{{stroke:var(--hair);stroke-width:1}}
.rmsring{{fill:none;stroke-width:1.2;stroke-dasharray:3 3;opacity:.9}}
.airy{{fill:none;stroke:var(--ink2);stroke-width:1}}
.scalebar{{stroke:var(--ink2);stroke-width:1.6}}
.scalelab{{font-family:var(--mono);font-size:8.5px;fill:var(--muted)}}
figcaption{{display:grid;grid-template-columns:1fr 1fr;gap:3px 12px;margin-top:9px;
  font-family:var(--mono);font-size:11px;font-variant-numeric:tabular-nums}}
figcaption div{{display:flex;justify-content:space-between;border-bottom:1px solid var(--hair);padding:2px 0}}
figcaption span{{color:var(--muted)}} figcaption b{{font-weight:600}}
.legend{{display:flex;gap:16px;flex-wrap:wrap;font-family:var(--mono);font-size:11px;color:var(--ink2);margin-top:14px}}
.legend i{{display:inline-block;width:20px;height:0;border-top:1.2px dashed currentColor;vertical-align:middle;margin-right:5px}}
.notes{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:14px;margin-top:8px}}
.note{{border-left:2px solid var(--accent);padding:4px 0 4px 14px}}
.note h3{{font-size:13.5px;margin:0 0 4px}} .note p{{margin:0;font-size:13px;color:var(--ink2)}}
.foot{{color:var(--muted);font-family:var(--mono);font-size:11px;margin-top:26px;
  border-top:1px solid var(--line);padding-top:14px}}
b.hl{{color:var(--ink)}}
</style>
<div class="wrap"><div class="inner">
  <div class="top">
    <div>
      <p class="eyebrow">JASPER · 500–1000 nm on TCD1304 · cylindrical field flattener</p>
      <h1>Spot diagram across the 500–1000&nbsp;nm band</h1>
      <p class="sub">BeamFour ray trace: 100&nbsp;µm fiber → 30&nbsp;µm slit → f=78&nbsp;mm collimator → 600&nbsp;l/mm grating →
        f=78&nbsp;mm camera → <b class="hl">cylindrical field flattener</b> → TCD1304 (29.1&nbsp;mm), detector tilted 6°.
        169 hexapolar pupil rays/wavelength at f/4.9. The cylinder corrects the dispersion-plane field curvature,
        flattening the response to ~2&nbsp;nm across the band. Airy disk shown as the diffraction floor.</p>
    </div>
    <button class="tog" id="tog" aria-label="Toggle colour theme">◐ theme</button>
  </div>

  <div class="tiles">
    <div class="tile"><div class="k">Plate scale</div><div class="v">{plate_scale(500):.0f}–{plate_scale(1000):.0f}<small> nm/mm (local)</small></div></div>
    <div class="tile"><div class="k">Spectral span</div><div class="v">{span:.2f}<small> mm</small></div></div>
    <div class="tile"><div class="k">Best Δλ</div><div class="v">{best["dlam"]:.1f}<small> nm @ {best["w"]}nm</small></div></div>
    <div class="tile"><div class="k">Airy radius</div><div class="v">~4<small> µm (floor)</small></div></div>
    <div class="tile"><div class="k">Detector fill</div><div class="v">{span/DET_MM*100:.0f}<small> %</small></div></div>
  </div>

  <div class="card">
    <h2>Where each wavelength lands</h2>
    <p class="cap">The grating sorts colour by position on the 14.336&nbsp;mm sensor — this strip is the spectrum itself.</p>
    {strip}
  </div>

  <div class="card">
    <h2>RMSx (dispersion axis) vs wavelength</h2>
    <p class="cap">Only the dispersion-axis blur sets spectral resolution, so RMSx is plotted, not the radial RMS.
      With the cylindrical flattener the response is <em>flat</em>: mid-band sits near best focus while 500 and
      1000&nbsp;nm carry the small residual the corrector and 6° tilt cannot fully remove.</p>
    {bar}
  </div>

  <div class="card">
    <h2>Spot-diagram matrix — common ±430 µm scale</h2>
    <p class="cap">Each panel: {minn}–{maxn} of 169 traced ray hits at the detector, referenced to the spot centroid.
      Dashed ring = RMS radius · centre dot = Airy disk (diffraction floor, ~4diffraction floor, ~2&nbsp;µmnbsp;µm — sub-pixel here).
      Δλ is the FWHM-based resolution element.</p>
    <div class="matrix">{panels}</div>
    <div class="legend">
      <span><i></i>RMS-radius ring</span>
      <span>● Airy disk (≈4 µm)</span>
      <span>chip = spot regime</span>
    </div>
  </div>

  <div class="card">
    <h2>What the spots tell you</h2>
    <div class="notes">
      <div class="note"><h3>The cylinder flattens the field</h3><p>The limiter was the camera's
        <em>dispersion-plane field curvature</em> — best focus swung ~2&nbsp;mm across the band. A plano-cylindrical
        N-BK7 flattener (f&nbsp;≈&nbsp;−70&nbsp;mm) near the detector corrects it, bringing the response to
        <b>~1.9&nbsp;nm mean / 2.9&nbsp;nm worst</b>.</p></div>
      <div class="note"><h3>Why cylindrical, not spherical</h3><p>Only the dispersion axis sets resolution, so the
        corrector needs power in that plane <em>only</em>. A cylinder flattens the tangential focal surface without
        touching the slit-height axis — avoiding the whole-spectrum magnification that makes a spherical flattener
        a wash on a nearly-full detector.</p></div>
      <div class="note"><h3>Band edges are the residual</h3><p>500 and 1000&nbsp;nm focus deepest, so they carry the
        last of the curvature (~2.9&nbsp;nm) while mid-band sits near 1&nbsp;nm. The 6° detector tilt balances the
        two ends; a single cylinder cannot remove the small chromatic remainder.</p></div>
      <div class="note"><h3>Geometry-limited</h3><p>RMSx is far larger than the ~4&nbsp;µm Airy disk at f/4.9, so
        diffraction is negligible; the blur is ray aberration plus the residual defocus the tilt cannot remove.</p></div>
      <div class="note"><h3>The path to ~1&nbsp;nm</h3><p>Each wavelength reaches ~0.5&nbsp;nm at <em>its own</em>
        focus, and best focus is a smooth parabola in wavelength. A 2-exposure focus bracket (detector at two
        positions ~1.7&nbsp;mm apart, stitched) reaches <b>~1&nbsp;nm</b> — demonstrated in simulation, pending build.</p></div>
      <div class="note"><h3>Not shown: diffraction</h3><p>BeamFour is geometric — the Airy circle is computed
        analytically (1.22·λ·f/#). For the true PSF at band centre you'd need a physical-optics tool.</p></div>
    </div>
  </div>

  <p class="foot">Source: BeamFour InOut trace · 500-1000nm cylinder design · 3042 rays (3 fields × 169 pupil × 6 λ) ·
    {surv} centre-field survivors · centre field shown · generated {__import__("datetime").date.today().isoformat()}</p>
</div></div>
<script>
(function(){{
  var r=document.documentElement,b=document.getElementById('tog');
  b.addEventListener('click',function(){{
    var cur=r.getAttribute('data-theme');
    if(!cur)cur=matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light';
    r.setAttribute('data-theme',cur==='dark'?'light':'dark');
  }});
}})();
</script>'''
open(os.path.join(DESIGN,"JASPER_spot_analysis.html"),"w",encoding="utf-8").write(html)
print("wrote JASPER_spot_analysis.html  (%d chars)"%len(html))
