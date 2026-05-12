"""Replot CS3 precision figure as Δ Precision@K vs k-core baseline.

The previous figure stacked 10 absolute Precision@K curves on top of the
k-core line; the gap was <5% so all curves looked identical at this zoom.
This script plots the *gain* Precision@K(1,s) - Precision@K(k-core) so the
y=0 line *is* the k-core baseline and every curve above it is a win.

Output: cs3_precision.pdf overwritten with the cleaner view, plus
        cs3_precision.png for inline previews.
"""
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["text.usetex"] = False
matplotlib.rcParams["font.family"] = "DejaVu Sans"
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent

# Restrict to the s values where the active set still covers K=10000 so the
# Precision@K comparison is apples-to-apples (s=15+ drops below 10k).
S_VALUES = [3, 5, 8, 10]
# K grid
KS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]

# Load communities and vertex coverage
communities = []
vertex_to_comm = defaultdict(set)
with open(ROOT / "com-dblp.top5000.cmty.txt") as f:
    for idx, line in enumerate(f):
        members = [int(x) for x in line.strip().split()]
        communities.append(members)
        for v in members:
            vertex_to_comm[v].add(idx)
covered = set(vertex_to_comm.keys())

# Load (1,s)-core values
cores = {}
for s in S_VALUES:
    df = pd.read_csv(ROOT / f"com-dblp-snap_def_s{s}.tsv", sep="\t", comment="#",
                     names=["vid", "core"], dtype={"vid": np.int64, "core": np.float64})
    cores[s] = df.set_index("vid")["core"].values

# Load k-core baseline (cached)
kc = np.load(ROOT / "kcore.npy") if (ROOT / "kcore.npy").exists() else None
if kc is None:
    raise SystemExit("kcore.npy missing; run analyze_cs3.py first to generate it")

def precision_at_K(top_set, covered):
    if not top_set: return 0.0
    return len(top_set & covered) / len(top_set)

# Compute baseline Precision@K curve once
prec_kc = {}
for K in KS:
    top = set(int(i) for i in np.argsort(kc)[-K:])
    prec_kc[K] = precision_at_K(top, covered)

# Compute (1,s) Precision@K and Δ vs baseline
delta = {}
prec_s = {}
for s in S_VALUES:
    cv = cores[s]
    delta[s] = {}
    prec_s[s] = {}
    for K in KS:
        nz = int((cv > 0).sum())
        k_eff = min(K, nz)
        if k_eff == 0:
            continue
        top = set(int(i) for i in np.argsort(cv)[-k_eff:])
        p = precision_at_K(top, covered)
        prec_s[s][K] = p
        delta[s][K] = p - prec_kc[K]

# Layout: 2-panel, left = absolute Precision (zoomed), right = gap over k-core
WIN  = "#1f4e9c"   # SPIN★ blue palette for ours
CND  = "#c0392b"   # red for baseline
PALETTE = ["#1f4e9c", "#3a8df0", "#1aa39a", "#a050c0"]  # 4 s values

fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2))

# Left: absolute precision, k-core dashed red, (1,s) coloured
ax = axes[0]
ax.plot(KS, [prec_kc[K] for K in KS], color=CND, marker="s", ms=6, lw=2.0,
        ls="--", mec=CND, mfc=CND, label=r"$k$-core (baseline)", zorder=2)
for s, c in zip(S_VALUES, PALETTE):
    xs = sorted(prec_s[s].keys())
    ys = [prec_s[s][K] for K in xs]
    ax.plot(xs, ys, color=c, marker="o", ms=5, lw=1.8,
            mec=c, mfc=c, label=fr"$(1,{s})$-core", zorder=3)
ax.set_xscale("log")
ax.set_xlabel(r"$K$ (top ranking)", fontsize=10)
ax.set_ylabel(r"Precision@$K$", fontsize=10)
ax.set_title("Absolute Precision@$K$", fontsize=10.5, fontweight="bold", color="#111")
ax.grid(True, which="major", color="#e8e8e8", lw=0.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(loc="lower right", frameon=False, fontsize=8)

# Right: GAP (1,s) - k-core, baseline = y=0, gap shaded above
ax = axes[1]
ax.axhline(0.0, color=CND, ls="--", lw=1.6, label=r"$k$-core baseline (gap $=0$)")
for s, c in zip(S_VALUES, PALETTE):
    xs = sorted(delta[s].keys())
    ys = [delta[s][K] for K in xs]
    ax.plot(xs, ys, color=c, marker="o", ms=6, lw=2.0,
            mec=c, mfc=c, label=fr"$(1,{s})$-core", zorder=3)
ax.set_xscale("log")
ax.set_xlabel(r"$K$ (top ranking)", fontsize=10)
ax.set_ylabel(r"$\Delta$ Precision@$K$ (vs $k$-core)", fontsize=10)
ax.set_title("Gain over $k$-core baseline", fontsize=10.5, fontweight="bold", color="#111")
ax.grid(True, which="major", color="#e8e8e8", lw=0.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(loc="upper left", frameon=False, fontsize=8)

fig.tight_layout()
out_pdf = ROOT / "cs3_precision.pdf"
out_png = ROOT / "cs3_precision.png"
fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02)
fig.savefig(out_png, dpi=180, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print(f"wrote {out_pdf}")
print(f"wrote {out_png}")

# Print a small summary so the prose can cite the gain
print("\nMax Δ Precision@K per s:")
for s in S_VALUES:
    if delta[s]:
        K_best = max(delta[s], key=delta[s].get)
        print(f"  (1,{s})-core: +{delta[s][K_best]*100:.2f} pp at K={K_best}")
