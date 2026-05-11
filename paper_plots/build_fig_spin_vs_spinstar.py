"""
Figure: SPIN vs SPIN* — direct iteration vs event-driven optimization.

Story: SPIN is the literal fixed-point solver; SPIN* is its event-driven
optimization. Both compute the same core values. SPIN* avoids redundant
full-graph passes by only revisiting vertices whose path counters changed.

The figure shows time per (graph, s) cell as paired horizontal bars
sorted by speedup factor (largest first). SPIN gray, SPIN* black.

Per skill: monochrome, top-row legend, no grid.
"""
import csv, math, statistics
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

ROOT     = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "paper_data"
OUT_DIR  = ROOT / "paper_plots"

# ---- SPIN (LocalH) times from bench_local_v4.csv (took_ms is peel-only-equiv) ----
spin = {}
with (DATA_DIR / "bench_local_v4.csv").open() as f:
    for r in csv.DictReader(f):
        if r["status"] != "OK": continue
        raw = r.get("took_ms") or r.get("wall_ms")
        try:
            spin[(r["graph"], int(r["s"]))] = float(raw)
        except (KeyError, ValueError, TypeError):
            continue

# ---- SPIN* times from main benchmark v3 (peel_ms) ----
spinstar = defaultdict(list)
with (DATA_DIR / "01_main_benchmark_v3.csv").open() as f:
    for r in csv.DictReader(f):
        if r["status"] != "OK": continue
        if r["algorithm"] != "Pure": continue
        t = r.get("peel_ms") or r.get("took_ms")
        if not t: continue
        try:
            spinstar[(r["graph"], int(r["s"]))].append(float(t))
        except ValueError:
            continue
spinstar = {k: statistics.median(v) for k, v in spinstar.items()}

# ---- pair, compute speedup, sort by speedup descending ----
all_pairs = []
for key, t_spin in spin.items():
    t_star = spinstar.get(key)
    if not t_star: continue
    all_pairs.append((key, t_spin, t_star, t_spin / t_star))
all_pairs.sort(key=lambda x: -x[3])

# Cap to top-K cells by speedup so the tornado stays readable.
# Spread across graphs: pick at most PER_GRAPH cells per graph so the
# tornado isn't dominated by one input. Aggregate stats are reported
# over the full set of paired cells.
TOP_K = 18
PER_GRAPH = 3
per_count = {}
pairs = []
for p in all_pairs:
    g = p[0][0]
    if per_count.get(g, 0) >= PER_GRAPH: continue
    per_count[g] = per_count.get(g, 0) + 1
    pairs.append(p)
    if len(pairs) >= TOP_K: break

GRAPH_DISPLAY = {
    "com-amazon.ungraph":      "com-amazon",
    "soc-pokec-relationships": "soc-pokec",
    "twitter_combined":        "twitter",
}

# ---- figure ----
WIN  = "#1f4e9c"   # SPIN★ blue (matches main exp figures)
DIM  = "#2c2c2c"   # SPIN near-black (matches main exp figures)
ANNO = "#c0392b"   # speedup label red

# figsize larger than rendered (figure* = 7.0"); LaTeX scales down ⇒
# text in figure renders slightly smaller than body text.
n = len(pairs)
fig = plt.figure(figsize=(9.0, 0.6 + 0.30 * n))
gs = GridSpec(2, 1, figure=fig, height_ratios=[0.10, 1.0],
              hspace=0.20, top=0.96, bottom=0.07, left=0.18, right=0.97)

# Top legend
ax_leg = fig.add_subplot(gs[0, 0]); ax_leg.axis("off")
handles = [
    Line2D([0], [0], color=WIN, lw=6, label="SPIN$^\\star$ (incremental)"),
    Line2D([0], [0], color=DIM, lw=6, label="SPIN (direct iteration)"),
]
ax_leg.legend(handles=handles, loc="center", ncol=2, frameon=False,
              fontsize=10, handlelength=2.2, columnspacing=2.5)

ax = fig.add_subplot(gs[1, 0])

ypos = list(range(len(pairs)))[::-1]
labels  = [f"{GRAPH_DISPLAY.get(g, g)} $s{{=}}{s}$" for (g, s), _, _, _ in pairs]
spin_t  = [p[1] for p in pairs]
star_t  = [p[2] for p in pairs]
speedup = [p[3] for p in pairs]

bar_h = 0.38
ax.barh([y + bar_h/2 for y in ypos], spin_t, bar_h, color=DIM,
        edgecolor="white", linewidth=0.7)
ax.barh([y - bar_h/2 for y in ypos], star_t, bar_h, color=WIN,
        edgecolor="white", linewidth=0.7)

# value labels at end of each bar
for y, ts, tt in zip(ypos, spin_t, star_t):
    ax.text(ts * 1.18, y + bar_h/2, f"{ts:,.1f} ms",
            va="center", ha="left", fontsize=8, color="#444")
    ax.text(tt * 1.18, y - bar_h/2, f"{tt:,.1f} ms",
            va="center", ha="left", fontsize=8, color="#222", fontweight="bold")

# speedup annotation in red, on the right of SPIN bar
xmax = max(spin_t) * 30
for y, sp in zip(ypos, speedup):
    if sp >= 100:
        txt = f"${sp:,.0f}\\times$"
    elif sp >= 10:
        txt = f"${sp:.1f}\\times$"
    else:
        txt = f"${sp:.2f}\\times$"
    ax.text(xmax * 0.55, y, txt,
            va="center", ha="left", fontsize=9.5,
            color=ANNO, fontweight="bold")

ax.set_yticks(ypos)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xscale("log")
ax.set_xlim(0.5, xmax)
ax.set_xlabel("time (ms, log scale)", fontsize=10)

# Header above the speedup column
ax.text(xmax * 0.55, ypos[0] + 1.0, "speedup",
        va="bottom", ha="left", fontsize=9.5, color=ANNO,
        fontweight="bold")

for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.tick_params(axis="x", labelsize=8.5, colors="#444")
ax.grid(False)

out_pdf = OUT_DIR / "fig_spin_vs_spinstar.pdf"
out_png = OUT_DIR / "fig_spin_vs_spinstar.png"
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"wrote {out_pdf}")
print(f"  shown (top-{TOP_K}): gmean = "
      f"{math.exp(sum(math.log(p[3]) for p in pairs) / len(pairs)):,.1f}x")
print(f"  full set ({len(all_pairs)} pairs): gmean = "
      f"{math.exp(sum(math.log(p[3]) for p in all_pairs) / len(all_pairs)):,.1f}x, "
      f"max = {max(p[3] for p in all_pairs):,.0f}x, "
      f"min = {min(p[3] for p in all_pairs):.2f}x")
