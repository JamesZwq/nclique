#!/usr/bin/env python3
"""
make_rn_scaling_fig.py — the size-free scaling figure: region-native
counting time must NOT grow significantly with r or s (the user's main
control criterion), in contrast to CND/RegND* whose time explodes with s.

Reads:
  paper_data/bench_region_native.csv  (region-native: rn_total_s per r,s)
  paper_data/bench_full_merged.csv    (CND=REF, RegND*=RegNDC total_ms)
Writes:
  Sigmod2027Nuclear/figures/exp_rn_scaling.pdf
and prints a terse flatness summary (per-graph min/max region-native time
across the whole (r,s) grid -> flatness ratio; CND range for contrast).
"""
import csv, os, math
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RN = os.path.join(ROOT, "paper_data", "bench_region_native.csv")
MG = os.path.join(ROOT, "paper_data", "bench_full_merged.csv")
OUT = os.path.join(ROOT, "Sigmod2027Nuclear", "figures", "exp_rn_scaling.pdf")
GRAPHS = ["dblp-core30", "ca-GrQc", "ca-HepPh", "ca-CondMat", "com-dblp", "web-it-2004"]

def f(x):
    try: return float(x)
    except: return None

# region-native: (graph,r,s) -> total_s
rn = {}
for x in csv.DictReader(open(RN)):
    t = f(x["rn_total_s"])
    if t is not None and t >= 0:
        rn[(x["graph"], int(x["r"]), int(x["s"]))] = t

# CND/RegND* totals (s in seconds) from merged
cnd, reg = {}, {}
if os.path.exists(MG):
    for x in csv.DictReader(open(MG)):
        if x["graph"] not in GRAPHS or x["status"] != "OK": continue
        t = f(x.get("total_ms"))
        if t is None or t < 0: continue
        k = (x["graph"], int(x["r"]), int(x["s"]))
        if x["algo"] == "REF": cnd[k] = t / 1000.0
        elif x["algo"] == "RegNDC": reg[k] = t / 1000.0

# terse flatness summary
print(f"{'graph':13} | region-native total_s over all (r,s)         | CND total_s")
print(f"{'':13} | {'min':>9}{'max':>11}{'flat ratio':>12}  cells | {'min':>8}{'max':>10}")
for g in GRAPHS:
    vs = [v for (gg, r, s), v in rn.items() if gg == g]
    cs = [v for (gg, r, s), v in cnd.items() if gg == g]
    if not vs:
        print(f"{g:13} | (no region-native data)"); continue
    fr = max(vs) / max(min(vs), 1e-6)
    cstr = f"{min(cs):8.2f}{max(cs):10.1f}" if cs else "     n/a"
    print(f"{g:13} | {min(vs):9.3f}{max(vs):11.2f}{fr:11.1f}x {len(vs):5} | {cstr}")

# figure: per-graph panel, x=s, region-native (r=3) vs CND (r=3) vs RegND* (r=3)
ncol = 3; nrow = 2
fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 3.1, nrow * 2.3))
R0 = 3
for i, g in enumerate(GRAPHS):
    ax = axes[i // ncol][i % ncol]
    def series(d):
        pts = sorted((s, v) for (gg, r, s), v in d.items() if gg == g and r == R0)
        return [s for s, _ in pts], [v for _, v in pts]
    sx, sy = series(rn); ax.plot(sx, sy, "o-", color="black", ms=3, label="region-native")
    cx, cy = series(cnd)
    if cx: ax.plot(cx, cy, "s--", color="gray", ms=3, mfc="none", label="CND")
    gx, gy = series(reg)
    if gx: ax.plot(gx, gy, "^:", color="gray", ms=3, mfc="none", label="RegND*")
    ax.set_yscale("log"); ax.set_title(g, fontsize=8, fontweight="bold")
    ax.set_xlabel("s (r=3)", fontsize=7); ax.tick_params(labelsize=6)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if i == 0: ax.set_ylabel("time (s)", fontsize=7)
handles, labels = axes[0][0].get_legend_handles_labels()
fig.legend(handles, labels, ncol=3, loc="upper center", fontsize=7, frameon=False)
fig.tight_layout(rect=[0, 0, 1, 0.95])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT)
print(f"\n[fig] {OUT}")
