#!/usr/bin/env python3
"""Plot per-phase time and memory breakdown from phase_breakdown.tsv.

Emits:
  figures/exp_breakdown_time.pdf       4-panel stacked time bars (4 graphs × s axis)
  figures/breakdown_memory_DBLP.pdf    stacked delta-RSS for com-dblp
  figures/breakdown_memory_web-it-2004.pdf  same for web-it-2004

Paper §exp-breakdown-time recognises 5 logical phases:
  Tree    = buildSDCT
  MCEnum  = MCEnum
  Support = Support
  Index   = Index
  Peel    = Peel
Other PhaseLogger phases (loadAndSort, preMutation, prepareGraph,
dispatch_total) are lumped into 'Other' for clarity.

§exp-breakdown-memory uses 4 incremental slices:
  Graph   = loadAndSort  delta_rss
  Tree    = buildSDCT    delta_rss
  Index   = MCEnum + Index delta_rss
  RClique = Support + Peel delta_rss
"""
from __future__ import annotations
import argparse, csv, sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
TSV_DEFAULT = ROOT / "paper_data" / "phase_breakdown.tsv"
OUT_DIR_DEFAULT = ROOT / "Sigmod2027Nuclear" / "figures"

# Paper alias for graphs
GLABEL = {
    "ca-GrQc": "GRQ", "ca-CondMat": "CON", "ca-HepPh": "HEPP",
    "dblp-core30": "DBC", "com-dblp": "DB", "web-it-2004": "WI",
}

TIME_PHASES = ["buildSDCT", "MCEnum", "Index", "Support", "Peel"]
TIME_LABELS = {
    "buildSDCT": "Tree", "MCEnum": "MCEnum", "Index": "Index",
    "Support": "Support", "Peel": "Peel",
}

# Monochrome fills with hatch differentiation
TIME_STYLE = {
    "Tree":    dict(color="0.85", hatch="//"),
    "MCEnum":  dict(color="0.65", hatch="\\\\"),
    "Index":   dict(color="0.45", hatch="xx"),
    "Support": dict(color="0.25", hatch=".."),
    "Peel":    dict(color="0.05", hatch=""),
}

MEM_GROUPS = {
    "Graph":   ["loadAndSort"],
    "Tree":    ["buildSDCT"],
    "Index":   ["MCEnum", "Index"],
    "RClique": ["Support", "Peel"],
}
MEM_STYLE = {
    "Graph":   dict(color="0.85", hatch="//"),
    "Tree":    dict(color="0.60", hatch="\\\\"),
    "Index":   dict(color="0.35", hatch="xx"),
    "RClique": dict(color="0.10", hatch=""),
}


def load_tsv(path):
    """Return cells[(graph, r, s)] = {phase: (duration_ms, rss_kb, delta_rss_kb)}."""
    cells = defaultdict(dict)
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                graph, rval, sval, *_ = row["meta"].split(",")
                key = (graph, int(rval), int(sval))
            except (ValueError, KeyError):
                continue
            cells[key][row["phase"]] = (
                float(row["duration_ms"]),
                int(row["rss_kb"]),
                int(row["delta_rss_kb"]),
            )
    return cells


# ---------------------------------------------------------------------------
def plot_time_breakdown(cells, out_path, graphs=None):
    """Stacked horizontal bars: one column of bars per graph, one bar per s."""
    if graphs is None:
        graphs = ["ca-HepPh", "com-dblp", "ca-CondMat", "web-it-2004"]
    ncols = len(graphs)
    fig, axes = plt.subplots(
        1, ncols, figsize=(2.3 * ncols, 2.4), squeeze=False,
        gridspec_kw=dict(wspace=0.30),
    )
    fig.subplots_adjust(top=0.78, bottom=0.18, left=0.08, right=0.99)
    axes = axes[0]

    for ax_idx, g in enumerate(graphs):
        ax = axes[ax_idx]
        cells_g = sorted(k for k in cells if k[0] == g)
        if not cells_g:
            ax.text(0.5, 0.5, "(no data)", transform=ax.transAxes,
                    ha="center", va="center")
            ax.set_title(GLABEL.get(g, g))
            continue
        s_vals = [s for _, _, s in cells_g]
        # 5 phases per cell, in seconds
        durs_ms = {p: [] for p in TIME_PHASES}
        for k in cells_g:
            for p in TIME_PHASES:
                durs_ms[p].append(cells[k].get(p, (0,0,0))[0])

        x = np.arange(len(s_vals))
        bottom = np.zeros(len(s_vals))
        for p in TIME_PHASES:
            paper_label = TIME_LABELS[p]
            heights = [d / 1000.0 for d in durs_ms[p]]
            ax.bar(x, heights, bottom=bottom, width=0.7,
                   label=paper_label,
                   **TIME_STYLE[paper_label],
                   edgecolor="black", linewidth=0.5)
            bottom += np.array(heights)

        ax.set_xticks(x)
        ax.set_xticklabels([f"$s{{=}}{s}$" for s in s_vals], fontsize=8)
        ax.set_title(GLABEL.get(g, g), fontsize=10)
        if ax_idx == 0:
            ax.set_ylabel("time (s)", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(direction="out", length=3, labelsize=8)

    # Legend strip on top
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=len(TIME_PHASES),
                   bbox_to_anchor=(0.5, 1.00), frameon=False, fontsize=9,
                   handlelength=2.4)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[time]   -> {out_path}")


# ---------------------------------------------------------------------------
def plot_memory_breakdown_one(cells, graph, out_path):
    """Stacked horizontal bar per s value: 4 slices Graph/Tree/Index/RClique."""
    cells_g = sorted(k for k in cells if k[0] == graph)
    if not cells_g:
        print(f"[mem]    skip {graph}: no data", file=sys.stderr)
        return

    s_vals = [s for _, _, s in cells_g]
    slice_mb = {grp: [] for grp in MEM_GROUPS}
    for k in cells_g:
        for grp, phases in MEM_GROUPS.items():
            tot = sum(max(0, cells[k].get(p, (0,0,0))[2]) for p in phases)
            slice_mb[grp].append(tot / 1024.0)

    fig, ax = plt.subplots(figsize=(2.5, 2.2))
    fig.subplots_adjust(left=0.20, bottom=0.28, top=0.95, right=0.97)

    x = np.arange(len(s_vals))
    bottom = np.zeros(len(s_vals))
    for grp in MEM_GROUPS:
        h = np.array(slice_mb[grp])
        ax.bar(x, h, bottom=bottom, width=0.7,
               label=grp, **MEM_STYLE[grp],
               edgecolor="black", linewidth=0.5)
        bottom += h

    ax.set_xticks(x)
    ax.set_xticklabels([f"$s{{=}}{s}$" for s in s_vals], fontsize=8)
    ax.set_ylabel("memory (MB)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, labelsize=8)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15),
              fontsize=7.5, frameon=False, ncol=4,
              handlelength=1.4, columnspacing=0.8)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[mem]    -> {out_path}  ({len(s_vals)} cells)")


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tsv", type=Path, default=TSV_DEFAULT)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR_DEFAULT)
    args = p.parse_args()

    if not args.tsv.exists():
        print(f"missing TSV: {args.tsv}", file=sys.stderr); sys.exit(1)

    cells = load_tsv(args.tsv)
    print(f"[load] {len(cells)} cells from {args.tsv}")

    plot_time_breakdown(cells, args.out_dir / "exp_breakdown_time.pdf")
    plot_memory_breakdown_one(cells, "com-dblp",
                              args.out_dir / "breakdown_memory_DBLP.pdf")
    plot_memory_breakdown_one(cells, "web-it-2004",
                              args.out_dir / "breakdown_memory_web-it-2004.pdf")


if __name__ == "__main__":
    main()
