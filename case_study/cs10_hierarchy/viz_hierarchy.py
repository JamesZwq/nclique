"""
CS10 viz: branch-persistence fingerprint per domain at s=5.

One panel, three lines (one per graph), log-log axes:
    x = branch rank (1..N), sorted by persistence descending
    y = persistence = k_birth - k_death of that branch

Reader interpretation:
    long flat tail  -> rich nested hierarchy  (e.g., com-dblp)
    sharp drop      -> shallow hierarchy      (e.g., com-youtube)

Data:
    cs10_hierarchy_tree.csv with columns
        graph,s,node_id,k_birth,k_death,size_birth,size_death,parent_id

Outputs:
    cs10_hierarchy.{pdf,png}
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"]  = 42
matplotlib.rcParams["text.usetex"]  = False
matplotlib.rcParams["font.family"]  = "serif"
import matplotlib.pyplot as plt

ROOT     = Path(__file__).parent
TREE_CSV = ROOT / "cs10_hierarchy_tree.csv"

GRAPHS = [
    ("com-dblp",     "#1f78b4", "o"),
    ("com-youtube",  "#e31a1c", "s"),
    ("web-Stanford", "#33a02c", "^"),
]
S_FOR_TREE = 5


def load_persistence():
    by_g = defaultdict(list)
    with open(TREE_CSV) as f:
        for r in csv.DictReader(f):
            if int(r["s"]) != S_FOR_TREE: continue
            p = float(r["k_birth"]) - float(r["k_death"])
            by_g[r["graph"]].append(p)
    for g in by_g:
        by_g[g].sort(reverse=True)
    return by_g


def main():
    data = load_persistence()

    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    for gname, color, marker in GRAPHS:
        ps = data.get(gname, [])
        if not ps: continue
        ranks  = list(range(1, len(ps) + 1))
        ys     = [max(p, 1) for p in ps]
        ax.plot(ranks, ys, color=color, linewidth=1.6,
                label=gname, alpha=0.9)
        sample_idx = [0, len(ranks)//8, len(ranks)//3,
                      len(ranks)//2, 3*len(ranks)//4]
        for idx in sample_idx:
            if 0 <= idx < len(ranks):
                ax.scatter([ranks[idx]], [ys[idx]],
                           s=22, c=color, marker=marker,
                           edgecolors="white", linewidths=0.5, zorder=3)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("branch rank  (sorted by persistence)", fontsize=10.5)
    ax.set_ylabel(r"persistence  $k_{\mathrm{birth}}-k_{\mathrm{death}}$",
                  fontsize=10.5)
    ax.tick_params(axis="both", which="major", labelsize=9)
    ax.grid(True, which="major", alpha=0.25, linestyle=":")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.legend(fontsize=10, frameon=False, loc="lower left")

    fig.tight_layout()
    pdf = ROOT / "cs10_hierarchy.pdf"
    png = ROOT / "cs10_hierarchy.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    print(f"wrote {pdf}\nwrote {png}")


if __name__ == "__main__":
    main()
