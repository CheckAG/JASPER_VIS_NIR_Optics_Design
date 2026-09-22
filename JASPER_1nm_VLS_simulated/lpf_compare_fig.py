#!/usr/bin/env python3
"""LPF resolution + astigmatism for the SIMULATED design, as a 3-step progression that shows what
each correction buys:
    1. Concave spherical, no VLS      -> poor resolution, poor throughput   (baseline)
    2. Concave toroidal, no VLS       -> poor resolution, GOOD throughput   (toroid fixes throughput)
    3. Concave toroidal + VLS         -> GOOD resolution, GOOD throughput    (VLS fixes resolution)
The point: resolution and throughput are decoupled — VLS fixes one, the toroid fixes the other.
Reads lpf_all_designs.json (run lpf_resolution.py first). Palette validated CVD-safe (ΔE 16)."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
D={d["name"]:d for d in json.load(open("lpf_all_designs.json"))}
sph=D["Simulated VLS — SPHERICAL"]; tor=D["Simulated VLS — TOROIDAL"]

# three configurations (label, resolution nm, astig mm, colour)
cfg=[
 ("Concave spherical\n(no VLS)",      sph["dlam_classical_mean"], sph["astig_mm"], "#e8912a"),
 ("Concave toroidal\n(no VLS)",       tor["dlam_classical_mean"], tor["astig_mm"], "#1aa5cf"),
 ("Concave toroidal\n+ VLS",          tor["dlam_corr_mean"],      tor["astig_mm"], "#2f8f5b"),
]
labels=[c[0] for c in cfg]; res=[c[1] for c in cfg]; ast=[c[2] for c in cfg]; cols=[c[3] for c in cfg]
INK="#222222"; MUTE="#5a6672"
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12.5,4.6))
y=np.arange(len(cfg))[::-1]

# ---- panel 1: resolution (linear) ----
ax1.barh(y,res,color=cols,height=0.62,zorder=3)
for yi,v in zip(y,res):
    ax1.text(v+0.08,yi,f"{v:.2f} nm",va="center",fontsize=10.5,fontweight="bold",color=INK)
ax1.set_yticks(y); ax1.set_yticklabels(labels,fontsize=9.5,color=INK)
ax1.set_xlim(0,4.4); ax1.set_ylim(-0.6,2.6)
ax1.set_xlabel("Spectral resolution Δλ (nm) — lower is better",fontsize=9.5,color=MUTE)
ax1.set_title("Resolution",fontsize=12,color=INK,fontweight="bold")
ax1.grid(axis="x",alpha=0.22); ax1.set_axisbelow(True)
for s in ("top","right","left"): ax1.spines[s].set_visible(False)
ax1.text(1.75,y[2],"◄  VLS fixes resolution",fontsize=9,color="#2f8f5b",va="center",ha="left",fontweight="bold")

# ---- panel 2: astigmatism / throughput (log) ----
ax2.barh(y,ast,color=cols,height=0.62,zorder=3)
for yi,a in zip(y,ast):
    cap=min(100,0.2/a*100) if a>0 else 100
    ax2.text(a*1.18,yi,f"{a:.2f} mm · ~{cap:.0f}% on pixel",va="center",fontsize=9.5,color=INK)
ax2.set_xscale("log"); ax2.set_yticks(y); ax2.set_yticklabels([]); ax2.set_xlim(0.02,60); ax2.set_ylim(-0.6,2.6)
ax2.axvline(0.2,ls="--",color=MUTE,lw=1); ax2.text(0.2,2.42,"200 µm pixel",fontsize=8,color=MUTE,ha="center")
ax2.set_xlabel("Astigmatic line height (mm) — throughput, not resolution",fontsize=9.5,color=MUTE)
ax2.set_title("Throughput",fontsize=12,color=INK,fontweight="bold")
ax2.grid(axis="x",which="both",alpha=0.18); ax2.set_axisbelow(True)
for s in ("top","right","left"): ax2.spines[s].set_visible(False)
ax2.text(1.2,1.0,"◄ toroid fixes throughput",fontsize=9,color="#1391bd",va="center",ha="left",fontweight="bold")

fig.suptitle("JASPER simulated VLS grating — what each correction buys (Noda light-path-function, f/5.3)",
             fontsize=12.5,y=1.03,color=INK)
fig.text(0.5,-0.06,"Resolution and throughput are independent knobs: the VLS groove correction sharpens the spectrum (3.7 → 0.9 nm), "
         "the toroidal substrate collects the light (1.0 → 0.04 mm astigmatic line, ~19% → ~100% on the pixel). "
         "Best-case correction; simulated / theoretical design.",fontsize=8,color=MUTE,ha="center",wrap=True)
fig.tight_layout(); fig.savefig("jasper_lpf_comparison.png",dpi=150,bbox_inches="tight")
print(f"wrote jasper_lpf_comparison.png | res {[round(r,2) for r in res]} | astig {[round(a,2) for a in ast]}")
