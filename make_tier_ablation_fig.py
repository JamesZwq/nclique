#!/usr/bin/env python3
"""
Tier ablation figure for the SIGMOD 2027 nuclear-decomposition paper.

Reads paper_data/bench_tier_ablation_results.csv and writes one PDF:
    figures/fig_tier_ablation.pdf

Layout (figure*, 2-column):
    one row per of the 6 paper graphs (PAPER6 ordering, sorted by |V|),
    one bar group per (r, s) cell inside the row,
    four colored bars per group for T1, T2, T3, T4 (wall-clock, log y).
    TIMEOUT  -> hatched bar at the 3600 s cap
    SKIP     -> 'x' glyph drawn at the bar base (T1 only)

Annotations:
    on top of each bar group, two stacked text labels report
        "T4/T3"  speedup ratio (if both tiers are OK)
        "T3/T2"  speedup ratio (if both tiers are OK)
    using the same fmt_speedup logic as scripts/make_sigmod_figs.py.

Legend:
    bottom strip explains tier semantics
        T1 = direct enumeration
        T2 = +CPI
        T3 = +Index Elimination (IE)
        T4 = +V-safe pruning
    plus TIMEOUT and SKIP_T1_HARD markers.

Style:
    Reuses the exact knobs from scripts/make_sigmod_figs.py:
        figsize width 7.05 in (figure*),
        dpi 230, savefig bbox_inches='tight',
        spines top/right hidden, no gridlines,
        compact_log_tick formatter on the y-axis,
        small monospaced/sans tick fonts in the 5.4..8.2 pt range,
        monochrome grayscale palette for the four tiers.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, NullFormatter


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
CSV_DEFAULT = ROOT / "paper_data" / "bench_tier_ablation_results.csv"
OUT_DEFAULT = Path(
    "/Users/zhangwenqian/Library/CloudStorage/"
    "Dropbox/应用/Overleaf/Sigmod2027Nuclear/figures/"
    "fig_tier_ablation.pdf"
)


# ---------------------------------------------------------------------------
# Constants -- mirror make_sigmod_figs.py
# ---------------------------------------------------------------------------
PAPER6 = [
    "dblp-core30",
    "ca-GrQc",
    "ca-HepPh",
    "ca-CondMat",
    "com-dblp",
    "web-it-2004",
]
GRAPH_SHORT = {
    "dblp-core30": "DBC",
    "ca-GrQc": "GRQ",
    "ca-HepPh": "HEPP",
    "ca-CondMat": "CON",
    "com-dblp": "DB",
    "web-it-2004": "WI",
}

TIER_ORDER = ["T1", "T2", "T3", "T4"]
TIER_LONG = {
    "T1": "T1: direct",
    "T2": "T2: +CPI",
    "T3": "T3: +IE",
    "T4": "T4: +V-safe",
}
# Monochrome grayscale palette (paper-architect rule: no rainbow).
# Darker = more pruning machinery.
TIER_COLOR = {
    "T1": "0.85",
    "T2": "0.55",
    "T3": "0.25",
    "T4": "0.05",
}
TIER_EDGE = "black"
TIMEOUT_CAP_S = 3600.0          # 1 h budget used by the bench harness
TIMEOUT_HATCH = "////"
TIMEOUT_FACE = "#404040"        # matches make_sigmod_figs.heatmap timeout color
SKIP_GLYPH = "×"          # multiplication sign (paper-safe "x" marker)


# ---------------------------------------------------------------------------
# Helpers borrowed from make_sigmod_figs.py
# ---------------------------------------------------------------------------
def compact_log_tick(x, _pos):
    if x <= 0:
        return ""
    if 0.1 <= x < 10:
        return f"{x:g}"
    exp = math.log10(x)
    if abs(exp - round(exp)) < 1e-9:
        return rf"$10^{{{int(round(exp))}}}$"
    return f"{x:g}"


def fmt_speedup(x: float) -> str:
    if not math.isfinite(x) or x <= 0:
        return ""
    if x >= 1000:
        return f"{x/1000:.1f}k$\\times$"
    if x >= 100:
        return f"{x:.0f}$\\times$"
    if x >= 10:
        return f"{x:.1f}$\\times$"
    return f"{x:.2f}$\\times$"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_rows(csv_path: Path) -> list[dict]:
    with open(csv_path) as f:
        return list(csv.DictReader(f))


def classify_status(raw: str) -> str:
    """Normalize the bench harness's status strings into 3 buckets."""
    s = (raw or "").strip()
    if s == "OK":
        return "OK"
    if s == "TIMEOUT(skip)":
        return "SKIP"
    if s.startswith("TIMEOUT"):
        return "TIMEOUT"
    # Anything else (OOM, ERROR) treated as TIMEOUT-style overlay so the
    # figure does not silently swallow the cell.
    return "TIMEOUT"


def index_rows(rows: list[dict]):
    """idx[(graph, r, s, tier)] = (status_bucket, wall_s_or_None)."""
    idx: dict[tuple[str, int, int, str], tuple[str, float | None]] = {}
    for row in rows:
        try:
            g = row["graph"]
            r = int(row["r"])
            s = int(row["s"])
            # CSV `tier` column is the int 1..4; the bar lookup uses "T1".."T4"
            # which the harness writes to the `algo` column.
            tier = row.get("algo") or f"T{row['tier']}"
        except (KeyError, ValueError):
            continue
        bucket = classify_status(row.get("status", ""))
        wall_s: float | None = None
        if bucket == "OK":
            try:
                ms = float(row["wall_ms"])
                if ms > 0:
                    wall_s = ms / 1000.0
            except (KeyError, ValueError, TypeError):
                wall_s = None
        idx[(g, r, s, tier)] = (bucket, wall_s)
    return idx


def collect_cells(idx) -> dict[str, list[tuple[int, int]]]:
    """Per-graph sorted list of distinct (r, s) cells present in the CSV."""
    per: dict[str, set[tuple[int, int]]] = defaultdict(set)
    for (g, r, s, _tier) in idx.keys():
        per[g].add((r, s))
    return {g: sorted(per[g]) for g in PAPER6 if per[g]}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def draw_cell(ax, idx, graph: str, r: int, s: int):
    """Render one (r, s) bar group for one graph onto a dedicated axes."""
    xs = np.arange(len(TIER_ORDER))
    width = 0.78

    bar_top: dict[str, float] = {}     # tier -> top (s) for annotation
    statuses: dict[str, str] = {}

    for xi, tier in enumerate(TIER_ORDER):
        bucket, wall_s = idx.get((graph, r, s, tier), ("MISSING", None))
        statuses[tier] = bucket

        if bucket == "OK" and wall_s is not None:
            ax.bar(
                xi, wall_s, width=width,
                color=TIER_COLOR[tier], edgecolor=TIER_EDGE, linewidth=0.45,
            )
            bar_top[tier] = wall_s
        elif bucket == "TIMEOUT":
            ax.bar(
                xi, TIMEOUT_CAP_S, width=width,
                facecolor=TIMEOUT_FACE, edgecolor=TIER_EDGE, linewidth=0.45,
                hatch=TIMEOUT_HATCH,
            )
            bar_top[tier] = TIMEOUT_CAP_S
        elif bucket == "SKIP":
            # Draw a low stub bar so the slot is visible, then mark with 'x'.
            ax.bar(
                xi, TIMEOUT_CAP_S, width=width,
                facecolor="white", edgecolor=TIER_EDGE, linewidth=0.45,
                hatch=TIMEOUT_HATCH,
            )
            ax.text(
                xi, TIMEOUT_CAP_S * 0.18, SKIP_GLYPH,
                ha="center", va="center", fontsize=8.2, color="black",
                fontweight="bold",
            )
            bar_top[tier] = TIMEOUT_CAP_S
        else:
            # Truly missing: leave the slot blank.
            bar_top[tier] = float("nan")

    # ---- speedup annotations (T4/T3 and T3/T2 when both OK) ----
    annotations: list[str] = []
    if statuses.get("T3") == "OK" and statuses.get("T4") == "OK":
        sp = bar_top["T3"] / bar_top["T4"] if bar_top["T4"] > 0 else float("nan")
        annotations.append(("T4/T3", fmt_speedup(sp)))
    if statuses.get("T2") == "OK" and statuses.get("T3") == "OK":
        sp = bar_top["T2"] / bar_top["T3"] if bar_top["T3"] > 0 else float("nan")
        annotations.append(("T3/T2", fmt_speedup(sp)))

    finite_tops = [v for v in bar_top.values() if math.isfinite(v) and v > 0]
    if finite_tops:
        anchor = max(finite_tops)
    else:
        anchor = TIMEOUT_CAP_S

    # Two-line annotation, stacked above the tallest bar.
    for k, (tag, txt) in enumerate(annotations):
        if not txt:
            continue
        y = anchor * (4.2 if k == 0 else 1.7)
        ax.text(
            1.5, y, f"{tag}={txt}",
            ha="center", va="bottom", fontsize=5.7, color="black",
        )

    # ---- axes cosmetics ----
    ax.set_yscale("log")
    ax.set_ylim(0.05, TIMEOUT_CAP_S * 18.0)
    ax.set_xticks(xs)
    ax.set_xticklabels(TIER_ORDER, fontsize=5.4)
    ax.yaxis.set_major_formatter(FuncFormatter(compact_log_tick))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="both", labelsize=5.4, length=2.0, pad=1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(f"$(r,s)=({r},{s})$", fontsize=6.6, pad=2.0)


def make_figure(csv_path: Path, out_path: Path) -> dict:
    rows = load_rows(csv_path)
    idx = index_rows(rows)
    cells = collect_cells(idx)

    # Determine per-row column count (each graph has the same # of cells in
    # the bench, but guard against asymmetry by padding with hidden axes).
    ncols = max((len(v) for v in cells.values()), default=0)
    nrows = len(PAPER6)
    if ncols == 0 or nrows == 0:
        raise RuntimeError("no tier-ablation cells parsed from CSV")

    # 2-column SIGMOD figure* width = 7.05 in (matches make_sigmod_figs.py).
    fig_w = 7.05
    fig_h = 1.32 * nrows + 0.55
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(fig_w, fig_h),
        squeeze=False,
        gridspec_kw=dict(wspace=0.36, hspace=0.78),
    )
    fig.subplots_adjust(left=0.075, right=0.995, top=0.93, bottom=0.10)

    rendered_cells = 0
    status_hist: Counter = Counter()

    for ri, graph in enumerate(PAPER6):
        graph_cells = cells.get(graph, [])
        for ci in range(ncols):
            ax = axes[ri][ci]
            if ci >= len(graph_cells):
                ax.set_visible(False)
                continue
            r, s = graph_cells[ci]
            draw_cell(ax, idx, graph, r, s)
            rendered_cells += 1
            for tier in TIER_ORDER:
                bucket, _ = idx.get((graph, r, s, tier), ("MISSING", None))
                status_hist[bucket] += 1

            if ci == 0:
                ax.set_ylabel(
                    f"{GRAPH_SHORT[graph]}\nwall (s)",
                    fontsize=6.9, labelpad=2.0,
                )

    # ---- legend strip (top, single row) ----
    tier_handles = [
        Patch(
            facecolor=TIER_COLOR[t], edgecolor=TIER_EDGE,
            linewidth=0.45, label=TIER_LONG[t],
        )
        for t in TIER_ORDER
    ]
    extra_handles = [
        Patch(
            facecolor=TIMEOUT_FACE, edgecolor=TIER_EDGE, linewidth=0.45,
            hatch=TIMEOUT_HATCH, label=f"TIMEOUT ({int(TIMEOUT_CAP_S)} s cap)",
        ),
        Line2D(
            [0], [0], marker="x", linestyle="None",
            markerfacecolor="black", markeredgecolor="black",
            markersize=7, label="SKIP_T1_HARD",
        ),
    ]
    fig.legend(
        handles=tier_handles + extra_handles,
        loc="upper center", ncol=6, frameon=False,
        fontsize=7.1, handlelength=1.6,
        bbox_to_anchor=(0.5, 1.00),
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=230, bbox_inches="tight")
    plt.close(fig)

    return {
        "rendered_cells": rendered_cells,
        "nrows": nrows,
        "ncols": ncols,
        "status_hist": dict(status_hist),
        "out_path": str(out_path),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path, default=CSV_DEFAULT)
    p.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = p.parse_args(argv)

    if not args.csv.exists():
        print(f"missing CSV: {args.csv}", file=sys.stderr)
        return 1

    info = make_figure(args.csv, args.out)
    # One-line debug summary -- spec requires this.
    hist = info["status_hist"]
    hist_str = ",".join(f"{k}={v}" for k, v in sorted(hist.items()))
    print(
        f"[tier-ablation] rendered {info['rendered_cells']} cells "
        f"(grid {info['nrows']}x{info['ncols']}) "
        f"statuses[{hist_str}] -> {info['out_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
