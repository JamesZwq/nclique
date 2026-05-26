#!/usr/bin/env python3
"""
Regenerate SIGMOD 2027 paper figures from merged bench CSV.

Reads paper_data/bench_v3_all_merged.csv and writes:
  - figures/exp_heatmap.pdf       Exp-1: time x coverage, 6 graphs x 3 algos
  - figures/exp_memory_r.pdf      Exp-2: M_REF / M_V3LM vs s, one curve per r

Paper aliases:
  REF    = Nucleus (Sariyuce et al.)
  RegNDC = Nuclear CD
  V3LM   = RegND (regnd) -- main algorithm

Paper-6 graphs (sorted by |V|):
  dblp-core30, ca-GrQc, ca-HepPh, ca-CondMat, com-dblp, web-it-2004

Plotting follows paper-architect rules: no gridlines, hidden top/right
spines, single legend row, monochrome where applicable.
"""
from __future__ import annotations
import argparse, csv, math, sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
CSV_DEFAULT = ROOT / "paper_data" / "bench_v3_all_merged.csv"
OUT_DEFAULT = ROOT / "Sigmod2027Nuclear" / "figures"

PAPER6 = ["dblp-core30", "ca-GrQc", "ca-HepPh", "ca-CondMat", "com-dblp", "web-it-2004"]
PAPER6_LABELS = {
    "dblp-core30": r"\textsc{dblp-core30}",
    "ca-GrQc":     r"\textsc{ca-GrQc}",
    "ca-HepPh":    r"\textsc{ca-HepPh}",
    "ca-CondMat":  r"\textsc{ca-CondMat}",
    "com-dblp":    r"\textsc{com-dblp}",
    "web-it-2004": r"\textsc{web-it-2004}",
}
# Plain labels (no LaTeX rendering at matplotlib level -- LaTeX compiles
# the caption, panel labels remain plain text).
PAPER6_PLAIN = {g: g for g in PAPER6}

ALGO_ALIAS = {"REF": "Nuclear CD", "RegNDC": "RegND"}
ALGOS_HEAT = ["REF", "RegNDC"]

# Paper's r-values for the memory ratio curves.
R_VALUES = [3, 4, 5, 7, 10, 15]

# Monochrome line styles (paper-architect rule: no rainbow).
LINE_STYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 2))]
MARKERS     = ["o", "s", "^", "D", "v", "x"]


# ---------------------------------------------------------------------------
def load_rows(csv_path: Path) -> list[dict]:
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def index_rows(rows: list[dict]) -> dict:
    """idx[(graph, r, s, algo)] = row."""
    idx = {}
    for row in rows:
        key = (row["graph"], int(row["r"]), int(row["s"]), row["algo"])
        idx[key] = row
    return idx


# ---------------------------------------------------------------------------
def make_heatmap(rows: list[dict], out_path: Path) -> None:
    """6 graphs x 3 algos heatmap. Cell color = log10(wall_ms / 1000).

    Special cells:
        TIMEOUT  -> dark gray
        OOM      -> patterned hatch
        N/A      -> white (no run attempted or r>=s)
    """
    idx = index_rows(rows)

    # r/s domain per graph (union over its REF/RegNDC/V3LM rows).
    domain = {}
    for g in PAPER6:
        rs = set()
        ss = set()
        for (gg, r, s, a), row in idx.items():
            if gg != g or a not in ALGOS_HEAT:
                continue
            rs.add(r); ss.add(s)
        if not rs:
            rs = {3}; ss = {4}
        domain[g] = (sorted(rs), sorted(ss))

    n_graphs = len(PAPER6)
    n_algos = len(ALGOS_HEAT)

    fig, axes = plt.subplots(
        n_algos, n_graphs,
        figsize=(2.3 * n_graphs, 1.7 * n_algos + 0.6),
        gridspec_kw=dict(wspace=0.20, hspace=0.32),
    )
    # Determine shared color scale (log seconds) across all OK cells.
    ok_seconds = []
    for row in rows:
        if row["status"] != "OK" or row["algo"] not in ALGOS_HEAT:
            continue
        if row["graph"] not in PAPER6:
            continue
        try:
            ok_seconds.append(float(row["wall_ms"]) / 1000.0)
        except (ValueError, KeyError):
            continue
    if not ok_seconds:
        print("[heatmap] no OK cells in CSV", file=sys.stderr); return

    vmin = max(min(ok_seconds), 1e-3)
    vmax = min(max(ok_seconds), 3600.0)
    norm = LogNorm(vmin=vmin, vmax=vmax)
    cmap = mpl.colormaps["YlOrRd"].copy()
    cmap.set_bad(color="white")

    for j, graph in enumerate(PAPER6):
        rs, ss = domain[graph]
        s_min, s_max = min(ss), max(ss)
        r_min, r_max = min(rs), max(rs)
        # Tight grid: row = r, col = s, only over [r_min..r_max] x [s_min..s_max]
        # but ensure r < s by masking the lower triangle.
        ax = axes[0, j] if n_algos > 1 else axes[j]
        for i, algo in enumerate(ALGOS_HEAT):
            ax = axes[i, j] if n_algos > 1 else axes[j]
            # Build matrix
            n_r = r_max - r_min + 1
            n_s = s_max - s_min + 1
            mat = np.full((n_r, n_s), np.nan)
            timeout_mask = np.zeros_like(mat, dtype=bool)
            oom_mask = np.zeros_like(mat, dtype=bool)
            for ri, r in enumerate(range(r_min, r_max + 1)):
                for si, s in enumerate(range(s_min, s_max + 1)):
                    if s <= r:
                        continue
                    row = idx.get((graph, r, s, algo))
                    if row is None:
                        continue
                    status = row["status"]
                    if status == "OK":
                        try:
                            mat[ri, si] = float(row["wall_ms"]) / 1000.0
                        except ValueError:
                            pass
                    elif status == "TIMEOUT":
                        timeout_mask[ri, si] = True
                    elif "OOM" in status:
                        oom_mask[ri, si] = True

            im = ax.imshow(
                mat,
                cmap=cmap, norm=norm,
                origin="lower",
                extent=(s_min - 0.5, s_max + 0.5, r_min - 0.5, r_max + 0.5),
                aspect="auto",
                interpolation="nearest",
            )
            # Overlay TIMEOUT (dark gray) and OOM (purple-ish gray)
            if timeout_mask.any():
                to_overlay = np.where(timeout_mask, 1.0, np.nan)
                ax.imshow(
                    to_overlay,
                    cmap=mpl.colors.ListedColormap(["#404040"]),
                    origin="lower",
                    extent=(s_min - 0.5, s_max + 0.5, r_min - 0.5, r_max + 0.5),
                    aspect="auto",
                    interpolation="nearest",
                )
            if oom_mask.any():
                oom_overlay = np.where(oom_mask, 1.0, np.nan)
                ax.imshow(
                    oom_overlay,
                    cmap=mpl.colors.ListedColormap(["#7a4f8a"]),
                    origin="lower",
                    extent=(s_min - 0.5, s_max + 0.5, r_min - 0.5, r_max + 0.5),
                    aspect="auto",
                    interpolation="nearest",
                )

            # Mask the r>=s region with light hatch to make it visually obvious
            for r in range(r_min, r_max + 1):
                ax.fill_between(
                    [s_min - 0.5, min(r + 0.5, s_max + 0.5)],
                    r - 0.5, r + 0.5,
                    color="#ececec", zorder=0,
                )

            # Cosmetics
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_color("#aaaaaa")
            ax.spines["bottom"].set_color("#aaaaaa")

            if j == 0:
                ax.set_ylabel(ALGO_ALIAS[algo], fontsize=10, rotation=0,
                              labelpad=22, ha="right", va="center")
            if i == 0:
                ax.set_title(PAPER6_PLAIN[graph], fontsize=10)
            if i == n_algos - 1:
                ax.set_xlabel("s", fontsize=9)
            # x-axis ticks: show first/last s
            ax.set_xticks([s_min, s_max])
            ax.set_xticklabels([str(s_min), str(s_max)], fontsize=7)
            if j == 0:
                ax.set_yticks([r_min, r_max])
                ax.set_yticklabels([str(r_min), str(r_max)], fontsize=7)

    # Colorbar
    cbar_ax = fig.add_axes([0.93, 0.18, 0.012, 0.68])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_label(r"time (s)", fontsize=9)
    cb.ax.tick_params(labelsize=7)

    # Legend strip for TIMEOUT / OOM / not-attempted
    legend_handles = [
        Patch(facecolor="#404040", label="timeout"),
        Patch(facecolor="#7a4f8a", label="OOM"),
        Patch(facecolor="#ececec", label=r"$r \geq s$"),
        Patch(facecolor="white", edgecolor="#aaaaaa", label="not attempted"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=4, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.45, -0.01))

    fig.subplots_adjust(left=0.10, right=0.91, top=0.92, bottom=0.13)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[heatmap] -> {out_path}")


# ---------------------------------------------------------------------------
def make_memory_ratio(rows: list[dict], out_path: Path) -> None:
    """Memory ratio M_REF / M_V3LM vs s, one curve per r in R_VALUES,
    one panel per paper-6 graph that has any matched cell."""
    idx = index_rows(rows)

    # Build per-graph data: data[graph][r] = list of (s, ratio).
    # Use RegNDC label (newer bench snapshot of the RegND algorithm);
    # fall back to V3LM for cells only present under the older label.
    data = defaultdict(lambda: defaultdict(list))
    for (g, r, s, a), row in idx.items():
        if g not in PAPER6 or a != "REF" or row["status"] != "OK":
            continue
        if r not in R_VALUES:
            continue
        v3 = idx.get((g, r, s, "RegNDC")) or idx.get((g, r, s, "V3LM"))
        if v3 is None or v3["status"] != "OK":
            continue
        try:
            m_ref = float(row["mem_kB"])
            m_v3 = float(v3["mem_kB"])
        except (ValueError, KeyError):
            continue
        if m_v3 <= 0:
            continue
        data[g][r].append((s, m_ref / m_v3))

    # Pick graphs that have at least one matched curve
    plotted_graphs = [g for g in PAPER6 if any(data[g][r] for r in R_VALUES)]
    if not plotted_graphs:
        print("[memory] no matched cells", file=sys.stderr); return

    n = len(plotted_graphs)
    cols = min(n, 4)
    rows_g = math.ceil(n / cols)
    fig, axes = plt.subplots(
        rows_g, cols,
        figsize=(2.4 * cols, 2.0 * rows_g + 0.4),
        squeeze=False,
        gridspec_kw=dict(wspace=0.32, hspace=0.42),
    )

    for idx_g, g in enumerate(plotted_graphs):
        ax = axes[idx_g // cols][idx_g % cols]
        for k, r in enumerate(R_VALUES):
            pts = sorted(data[g][r])
            if not pts:
                continue
            xs = [s for s, _ in pts]
            ys = [v for _, v in pts]
            ax.plot(xs, ys,
                    color="black",
                    linestyle=LINE_STYLES[k % len(LINE_STYLES)],
                    marker=MARKERS[k % len(MARKERS)],
                    markersize=4, linewidth=1.1,
                    markerfacecolor="white",
                    markeredgewidth=0.9,
                    label=f"r={r}")
        ax.axhline(1.0, color="#888888", linewidth=0.6, linestyle=":")
        ax.set_yscale("log")
        ax.set_title(PAPER6_PLAIN[g], fontsize=10)
        ax.set_xlabel("s", fontsize=9)
        if idx_g % cols == 0:
            ax.set_ylabel(r"$M_{\mathrm{NuclearCD}} / M_{\mathrm{RegND}}$", fontsize=9)
        ax.tick_params(axis="both", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Hide unused panels
    for k in range(len(plotted_graphs), rows_g * cols):
        axes[k // cols][k % cols].set_visible(False)

    # Single-row legend across the top
    handles = [
        Line2D([0], [0], color="black",
               linestyle=LINE_STYLES[k % len(LINE_STYLES)],
               marker=MARKERS[k % len(MARKERS)],
               markersize=4, linewidth=1.1,
               markerfacecolor="white", markeredgewidth=0.9,
               label=f"r={r}")
        for k, r in enumerate(R_VALUES)
    ]
    fig.legend(handles=handles, loc="upper center", ncol=len(R_VALUES),
               fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 1.02))

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[memory ] -> {out_path}")


# ---------------------------------------------------------------------------
def numerical_summary(rows: list[dict]) -> None:
    """Print the key numbers the paper text claims, recomputed from data."""
    idx = index_rows(rows)
    matched = []
    for (g, r, s, a), row in idx.items():
        if g not in PAPER6 or a != "REF" or row["status"] != "OK":
            continue
        v3 = idx.get((g, r, s, "RegNDC")) or idx.get((g, r, s, "V3LM"))
        if v3 is None or v3["status"] != "OK":
            continue
        try:
            t_ref = float(row["wall_ms"])
            t_v3  = float(v3["wall_ms"])
        except (ValueError, KeyError):
            continue
        if t_v3 <= 0 or t_ref <= 0:
            continue
        matched.append((g, r, s, t_ref / t_v3))

    print("\n=== Recomputed paper numbers ===")
    print(f"  matched (REF OK & RegND OK) cells: {len(matched)}")
    def gmean(xs):
        xs = [x for x in xs if x > 0]
        if not xs: return float("nan")
        return math.exp(sum(math.log(x) for x in xs) / len(xs))

    print(f"  geomean speedup (all 6 graphs): {gmean([x for *_, x in matched]):.2f}x")
    non_hepp = [x for g, *_ , x in matched if g != "ca-HepPh"]
    print(f"  geomean speedup (non-HepPh):    {gmean(non_hepp):.2f}x")
    dense = [x for g, *_ , x in matched if g in ("dblp-core30", "com-dblp", "web-it-2004")]
    print(f"  geomean speedup (dense 3):       {gmean(dense):.2f}x")
    dense_high_s = [x for g, r, s, x in matched
                    if g in ("dblp-core30", "com-dblp", "web-it-2004") and s >= 6]
    print(f"  geomean speedup (dense 3, s>=6): {gmean(dense_high_s):.2f}x")

    if matched:
        top = sorted(matched, key=lambda t: -t[3])[:5]
        print("  top-5 per-cell speedups:")
        for g, r, s, x in top:
            print(f"    {g:13s} ({r:2d},{s:2d}): {x:.1f}x")


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", type=Path, default=CSV_DEFAULT)
    p.add_argument("--out-dir", type=Path, default=OUT_DEFAULT)
    p.add_argument("--only", nargs="*", choices=["heatmap", "memory"],
                   default=["heatmap", "memory"])
    p.add_argument("--no-summary", action="store_true")
    args = p.parse_args()

    if not args.csv.exists():
        print(f"missing CSV: {args.csv}", file=sys.stderr); sys.exit(1)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.csv)
    print(f"[load] {len(rows)} rows from {args.csv}")

    if "heatmap" in args.only:
        make_heatmap(rows, args.out_dir / "exp_heatmap.pdf")
    if "memory" in args.only:
        make_memory_ratio(rows, args.out_dir / "exp_memory_r.pdf")
    if not args.no_summary:
        numerical_summary(rows)


if __name__ == "__main__":
    main()
