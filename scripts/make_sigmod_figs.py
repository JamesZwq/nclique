#!/usr/bin/env python3
"""
Regenerate SIGMOD 2027 paper figures from merged bench CSV.

Reads paper_data/bench_full_merged.csv and writes:
  - figures/exp_main_results.pdf   headline runtime/memory/coverage summary
  - figures/exp_time_advantage.pdf runtime completion profile
  - figures/exp_memory_advantage.pdf memory advantage by graph/setting
  - figures/exp_runtime.pdf       runtime on hard matched cells
  - figures/exp_memory_compare.pdf memory on the same hard cells
  - figures/exp_coverage.pdf      completed-cell coverage
  - figures/exp_scale_regnd.pdf   RegND-only scaling with s and r
  - figures/exp_heatmap.pdf       legacy full-grid time x coverage heatmap
  - figures/exp_memory_r.pdf      legacy memory-ratio sweep

Paper aliases:
  REF    = CND
  RegNDC = RegND*

Paper-6 graphs (sorted by |V|):
  dblp-core30, ca-GrQc, ca-HepPh, ca-CondMat, com-dblp, web-it-2004

Plotting follows paper-architect rules: no gridlines, hidden top/right
spines, single legend row, monochrome where applicable.
"""
from __future__ import annotations
import argparse, csv, math, sys
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullFormatter


# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
CSV_DEFAULT = ROOT / "paper_data" / "bench_full_merged.csv"
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
GRAPH_SHORT = {
    "dblp-core30": "DBC",
    "ca-GrQc": "GRQ",
    "ca-HepPh": "HEPP",
    "ca-CondMat": "CON",
    "com-dblp": "DB",
    "web-it-2004": "WI",
}

ALGO_ALIAS = {"REF": "CND", "RegNDC": "RegND*"}
ALGOS_HEAT = ["REF", "RegNDC"]
FORMAL_ALGOS = [
    ("REF", "CND", "0.50", "--", "s"),
    ("RegNDC", "RegND*", "black", "-", "o"),
    ("V3LM_NOCPI", "RegND", "0.66", ":", "^"),
    ("V3LM_HIER", "RegND-H", "0.22", "-.", "D"),
]
ALGO_CODE_ALIASES = {
    "REF": ("REF",),
    "RegNDC": ("RegNDC", "V3LM"),
    "V3LM_NOCPI": ("V3LM_NOCPI",),
    "V3LM_HIER": ("V3LM_HIER",),
}

# Paper's r-values for the memory ratio curves.
R_VALUES = [3, 4, 5, 7, 10, 15]
HARD_CASES = [
    ("dblp-core30", 5, 40),
    ("ca-GrQc", 7, 9),
    ("ca-CondMat", 7, 21),
    ("com-dblp", 5, 13),
]
DENSE_COVERAGE_GRAPHS = ["dblp-core30", "com-dblp", "web-it-2004"]
ADVANTAGE_R_ROWS = [4, 5, 6, 7]
ADVANTAGE_GRID_GRAPHS = [
    "dblp-core30",
    "ca-GrQc",
    "ca-CondMat",
    "com-dblp",
    "web-it-2004",
]

# Monochrome line styles (paper-architect rule: no rainbow).
LINE_STYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 2))]
MARKERS     = ["o", "s", "^", "D", "v", "x"]


def compact_log_tick(x, _pos):
    if x <= 0:
        return ""
    if 0.1 <= x < 10:
        return f"{x:g}"
    exp = math.log10(x)
    if abs(exp - round(exp)) < 1e-9:
        return rf"$10^{{{int(round(exp))}}}$"
    return f"{x:g}"


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


def regnd_row(idx: dict, graph: str, r: int, s: int):
    return idx.get((graph, r, s, "RegNDC")) or idx.get((graph, r, s, "V3LM"))


def row_for_algo(idx: dict, graph: str, r: int, s: int, algo: str):
    for code in ALGO_CODE_ALIASES.get(algo, (algo,)):
        row = idx.get((graph, r, s, code))
        if row is not None:
            return row
    return None


def to_float(row: dict, key: str):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return None


def ok_time_seconds(row: dict):
    if row is None or row.get("status") != "OK":
        return None
    ms = to_float(row, "wall_ms")
    return None if ms is None or ms <= 0 else ms / 1000.0


def ok_memory_gb(row: dict):
    if row is None or row.get("status") != "OK":
        return None
    kb = to_float(row, "mem_kB")
    return None if kb is None or kb <= 0 else kb / (1024.0 * 1024.0)


def fmt_speedup(x: float) -> str:
    if x >= 1000:
        return f"{x/1000:.1f}k$\\times$"
    if x >= 100:
        return f"{x:.0f}$\\times$"
    return f"{x:.1f}$\\times$"


def gmean(xs: list[float]) -> float:
    xs = [x for x in xs if x > 0 and math.isfinite(x)]
    if not xs:
        return float("nan")
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def collect_matched_metrics(rows: list[dict]) -> list[dict]:
    """All cells where both CND and RegND finish."""
    idx = index_rows(rows)
    matched = []
    for (graph, r, s, algo), ref in idx.items():
        if graph not in PAPER6 or algo != "REF" or ref.get("status") != "OK":
            continue
        reg = regnd_row(idx, graph, r, s)
        if reg is None or reg.get("status") != "OK":
            continue
        t_ref = to_float(ref, "wall_ms")
        t_reg = to_float(reg, "wall_ms")
        m_ref = to_float(ref, "mem_kB")
        m_reg = to_float(reg, "mem_kB")
        if not t_ref or not t_reg or not m_ref or not m_reg:
            continue
        matched.append(dict(
            graph=graph, r=r, s=s,
            time_ratio=t_ref / t_reg,
            memory_ratio=m_ref / m_reg,
        ))
    return matched


def add_count_labels(ax, xs, ys, counts, log=True):
    for x, y, n in zip(xs, ys, counts):
        if not math.isfinite(y) or y <= 0:
            continue
        if log:
            ax.text(x, y * 1.22, f"n={n}", ha="center", va="bottom",
                    fontsize=5.7)
        else:
            ax.text(x, y + 0.8, f"n={n}", ha="center", va="bottom",
                    fontsize=5.7)


def metric_value(row: dict | None, metric: str):
    if metric == "time":
        return ok_time_seconds(row)
    if metric == "memory":
        return ok_memory_gb(row)
    raise ValueError(metric)


def attempted_series_for(idx: dict, graph: str, r: int, algo: str, metric: str):
    codes = set(ALGO_CODE_ALIASES.get(algo, (algo,)))
    s_values = sorted({
        s for (g, rr, s, aa), row in idx.items()
        if g == graph and rr == r and aa in codes
    })
    rows_for_s = [(s, row_for_algo(idx, graph, r, s, algo)) for s in s_values]
    return [(s, metric_value(row, metric)) for s, row in rows_for_s]


def series_xy_with_gaps(series: list[tuple[int, float | None]]):
    """Break lines at failed rows and at missing s-values."""
    xs, ys = [], []
    prev_s = None
    for s, value in series:
        if prev_s is not None and s != prev_s + 1:
            xs.append(np.nan)
            ys.append(np.nan)
        if value is None:
            xs.append(s)
            ys.append(np.nan)
        else:
            xs.append(s)
            ys.append(value)
        prev_s = s
    return xs, ys


def completed_points_for(idx: dict, graph: str, r: int, algo: str, metric: str):
    pts = []
    for s, val in attempted_series_for(idx, graph, r, algo, metric):
        if val is not None:
            pts.append((s, val))
    return pts


def x_ticks_for(max_s: int):
    if max_s <= 35:
        ticks = [5, 15, max_s]
    elif max_s <= 55:
        ticks = [5, 20, max_s]
    elif max_s <= 130:
        ticks = [5, 40, 80, max_s]
    else:
        ticks = [5, 100, 200, 300, max_s]
    out = []
    for t in ticks:
        if t not in out:
            out.append(t)
    return out


def make_algorithm_grid(rows: list[dict], out_path: Path, metric: str, ylabel: str, tag: str) -> None:
    idx = index_rows(rows)
    nrows = len(ADVANTAGE_R_ROWS)
    ncols = len(ADVANTAGE_GRID_GRAPHS)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(7.05, 5.65),
        gridspec_kw=dict(wspace=0.25, hspace=0.34),
    )
    fig.subplots_adjust(left=0.065, right=0.995, top=0.905, bottom=0.075)

    for col, graph in enumerate(ADVANTAGE_GRID_GRAPHS):
        axes[0, col].set_title(GRAPH_SHORT[graph], fontsize=8.2, pad=3)

    graph_max_s = {}
    for graph in ADVANTAGE_GRID_GRAPHS:
        vals = []
        for r in ADVANTAGE_R_ROWS:
            for algo, _label, _color, _linestyle, _marker in FORMAL_ALGOS:
                vals.extend([s for s, _ in completed_points_for(idx, graph, r, algo, metric)])
        graph_max_s[graph] = max(vals) if vals else 10

    for row_i, r in enumerate(ADVANTAGE_R_ROWS):
        for col, graph in enumerate(ADVANTAGE_GRID_GRAPHS):
            ax = axes[row_i, col]
            any_pts = False
            for algo, label, color, linestyle, marker in FORMAL_ALGOS:
                series = attempted_series_for(idx, graph, r, algo, metric)
                pts = [(s, v) for s, v in series if v is not None]
                if not pts:
                    continue
                any_pts = True
                xs, ys = series_xy_with_gaps(series)
                ax.plot(
                    xs, ys,
                    color=color, linestyle=linestyle, linewidth=0.95,
                    marker=marker, markersize=1.9,
                    markevery=max(1, len(series) // 8),
                    markerfacecolor="white", markeredgewidth=0.5,
                    label=label,
                )
            if not any_pts:
                ax.text(0.5, 0.5, "no data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=6.2, color="0.45")

            max_s = graph_max_s[graph]
            ticks = x_ticks_for(max_s)
            ax.set_xlim(4, max_s * 1.03)
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(t) for t in ticks])
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(FuncFormatter(compact_log_tick))
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=5.4, length=2.0, pad=1.0)

            if col == 0:
                ax.set_ylabel(f"$r={r}$\n{ylabel}", fontsize=6.9)
            if row_i == nrows - 1:
                ax.set_xlabel("$s$", fontsize=6.9)
            else:
                ax.tick_params(labelbottom=False)

    handles = [
        Line2D([0], [0], color=color, linestyle=linestyle, marker=marker,
               markerfacecolor="white", markersize=3.0, linewidth=1.0,
               label=label)
        for _algo, label, color, linestyle, marker in FORMAL_ALGOS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               fontsize=7.1, bbox_to_anchor=(0.52, 0.995))
    fig.savefig(out_path, dpi=230, bbox_inches="tight")
    plt.close(fig)
    print(f"[{tag:<8}] -> {out_path}")


def make_time_advantage(rows: list[dict], out_path: Path) -> None:
    """Sorted per-configuration runtime (cactus plot) at fixed r.

    Cactus-plot convention (SAT/solver benchmarking): the y-axis is
    per-configuration runtime on a log scale, so a LOWER curve is faster
    and a curve that extends FARTHER RIGHT finishes more configurations.
    The lower-and-longer the curve, the better, which matches the
    "lower is better" reading readers expect from a runtime figure.
    This is the axis transpose of a completion profile.
    """
    idx = index_rows(rows)
    ncols = len(ADVANTAGE_R_ROWS)
    fig, axes = plt.subplots(
        1, ncols, figsize=(7.05, 2.15),
        gridspec_kw=dict(wspace=0.32),
    )
    fig.subplots_adjust(left=0.07, right=0.995, top=0.72, bottom=0.24)

    ymin, ymax = float("inf"), 0.0
    for ax, r in zip(axes, ADVANTAGE_R_ROWS):
        formal_codes = {
            code
            for algo, _label, _color, _linestyle, _marker in FORMAL_ALGOS
            for code in ALGO_CODE_ALIASES.get(algo, (algo,))
        }
        configs = sorted({
            (row["graph"], int(row["s"]))
            for row in rows
            if row["graph"] in ADVANTAGE_GRID_GRAPHS
            and int(row["r"]) == r
            and row["algo"] in formal_codes
        })
        denom = len(configs)
        for algo, label, color, linestyle, marker in FORMAL_ALGOS:
            times = []
            for graph, s in configs:
                row = row_for_algo(idx, graph, r, s, algo)
                t = ok_time_seconds(row)
                if t is not None:
                    times.append(t)
            times.sort()
            if not times:
                continue
            # x = cumulative fraction of configs solved (%), y = sorted
            # runtime (s). Each config is solved in non-decreasing time, so
            # the curve is monotone increasing left-to-right; a lower curve
            # solves the same fraction faster.
            xs = [100.0 * (i + 1) / denom for i in range(len(times))]
            ymin = min(ymin, times[0])
            ymax = max(ymax, times[-1])
            ax.step(xs, times, where="pre", color=color, linestyle=linestyle,
                    linewidth=1.2, label=label)
            step = max(1, len(times) // 8)
            ax.plot(
                xs[::step], times[::step], linestyle="None",
                marker=marker, color=color, markersize=2.4,
                markerfacecolor="white", markeredgewidth=0.6,
            )

        ax.set_yscale("log")
        ax.set_xlim(0, 103)
        ax.set_title(f"$r={r}$", fontsize=8.8, pad=3)
        ax.set_xlabel("configs solved (%)", fontsize=7.3)
        ax.tick_params(labelsize=6.2, length=2.0, pad=1.0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax is axes[0]:
            ax.set_ylabel("runtime (s)", fontsize=7.5)

    # Shared y-limits so the four panels are comparable; pad one decade.
    if ymax > 0 and math.isfinite(ymin):
        lo = 10 ** math.floor(math.log10(max(ymin, 1e-3)))
        hi = 10 ** math.ceil(math.log10(ymax))
        for ax in axes:
            ax.set_ylim(lo, hi)

    handles = [
        Line2D([0], [0], color=color, linestyle=linestyle, marker=marker,
               markerfacecolor="white", markersize=3.0, linewidth=1.1,
               label=label)
        for _algo, label, color, linestyle, marker in FORMAL_ALGOS
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               fontsize=7.0, bbox_to_anchor=(0.52, 0.99))
    fig.savefig(out_path, dpi=230, bbox_inches="tight")
    plt.close(fig)
    print(f"[timeprof] -> {out_path}")


def make_memory_advantage(rows: list[dict], out_path: Path) -> None:
    """Memory curves by fixed r and graph."""
    make_algorithm_grid(rows, out_path, "memory", "memory (GB)", "memadv")


# ---------------------------------------------------------------------------
def collect_hard_cases(rows: list[dict]) -> dict:
    idx = index_rows(rows)
    labels = []
    ncd_time = []
    reg_time = []
    ncd_memory = []
    reg_memory = []
    speedups = []
    memory_ratios = []

    for graph, r, s in HARD_CASES:
        ref = idx.get((graph, r, s, "REF"))
        reg = regnd_row(idx, graph, r, s)
        t_ref = ok_time_seconds(ref)
        t_reg = ok_time_seconds(reg)
        m_ref = ok_memory_gb(ref)
        m_reg = ok_memory_gb(reg)
        if t_ref is None or t_reg is None or m_ref is None or m_reg is None:
            raise RuntimeError(f"missing hard comparison cell: {graph} ({r},{s})")
        labels.append(f"{GRAPH_SHORT[graph]}\n({r},{s})")
        ncd_time.append(t_ref)
        reg_time.append(t_reg)
        ncd_memory.append(m_ref)
        reg_memory.append(m_reg)
        speedups.append(t_ref / t_reg)
        memory_ratios.append(m_ref / m_reg)

    return dict(
        labels=labels,
        ncd_time=ncd_time,
        reg_time=reg_time,
        ncd_memory=ncd_memory,
        reg_memory=reg_memory,
        speedups=speedups,
        memory_ratios=memory_ratios,
    )


def make_main_results(rows: list[dict], out_path: Path) -> None:
    """One compact headline figure: speed, memory, and completed cells."""
    hard = collect_hard_cases(rows)

    coverage_labels = [GRAPH_SHORT[g] for g in DENSE_COVERAGE_GRAPHS]
    coverage = {"CND": [], "RegND*": []}
    coverage_text = {"CND": [], "RegND*": []}
    for graph in DENSE_COVERAGE_GRAPHS:
        for algo, label in [("REF", "CND"), ("RegNDC", "RegND*")]:
            rows_ga = [row for row in rows if row["graph"] == graph and row["algo"] == algo]
            ok = sum(1 for row in rows_ga if row["status"] == "OK")
            total = len(rows_ga)
            pct = 100.0 * ok / total if total else 0.0
            coverage[label].append(pct)
            coverage_text[label].append(f"{ok:,}/{total:,}")

    fig, axes = plt.subplots(
        1, 3, figsize=(7.05, 1.95),
        gridspec_kw=dict(width_ratios=[1.18, 1.18, 1.0], wspace=0.40),
    )
    fig.subplots_adjust(left=0.07, right=0.995, top=0.77, bottom=0.28)

    cnd_style = dict(color="0.76", hatch="//", edgecolor="black", linewidth=0.5)
    reg_style = dict(color="0.08", edgecolor="black", linewidth=0.5)
    width = 0.34

    def grouped_bars(ax, labels, left_vals, right_vals, ylabel, title, annotations, log=True):
        x = np.arange(len(labels))
        ax.bar(x - width / 2, left_vals, width=width, label="CND", **cnd_style)
        ax.bar(x + width / 2, right_vals, width=width, label="RegND*", **reg_style)
        if log:
            ax.set_yscale("log")
            ymin = min(v for v in left_vals + right_vals if v > 0)
            ymax = max(left_vals + right_vals)
            ax.set_ylim(ymin * 0.55, ymax * 6.0)
        ax.set_title(title, fontsize=8.8, pad=3)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6.9)
        ax.tick_params(labelsize=7, length=2.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for i, txt in enumerate(annotations):
            y = max(left_vals[i], right_vals[i])
            ax.text(i, y * (1.28 if log else 1.03), txt,
                    ha="center", va="bottom", fontsize=6.3)

    grouped_bars(
        axes[0], hard["labels"], hard["ncd_time"], hard["reg_time"],
        "runtime (s)", "(a) time", [fmt_speedup(x) for x in hard["speedups"]],
    )
    grouped_bars(
        axes[1], hard["labels"], hard["ncd_memory"], hard["reg_memory"],
        "memory (GB)", "(b) memory", [fmt_speedup(x) for x in hard["memory_ratios"]],
    )

    ax = axes[2]
    x = np.arange(len(coverage_labels))
    ax.bar(x - width / 2, coverage["CND"], width=width, label="CND", **cnd_style)
    ax.bar(x + width / 2, coverage["RegND*"], width=width, label="RegND*", **reg_style)
    ax.set_ylim(0, 108)
    ax.set_title("(c) coverage", fontsize=8.8, pad=3)
    ax.set_ylabel("completed cells (%)", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(coverage_labels, fontsize=7)
    ax.tick_params(labelsize=7, length=2.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for i in range(len(coverage_labels)):
        ax.text(x[i] - width / 2, max(coverage["CND"][i] + 4, 5),
                coverage_text["CND"][i], ha="center", va="bottom",
                rotation=90, fontsize=5.3)
        ax.text(x[i] + width / 2, min(coverage["RegND*"][i] - 3, 96),
                coverage_text["RegND*"][i], ha="center", va="top",
                rotation=90, fontsize=5.3, color="white")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2,
               frameon=False, fontsize=8.3, handlelength=2.0,
               bbox_to_anchor=(0.52, 0.99))
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[main    ] -> {out_path}")


def make_runtime_compare(rows: list[dict], out_path: Path) -> None:
    """Runtime comparison on the strongest matched speedup cells."""
    data = collect_hard_cases(rows)
    labels = data["labels"]
    ncd_time = data["ncd_time"]
    reg_time = data["reg_time"]
    speedups = data["speedups"]

    fig, ax = plt.subplots(1, 1, figsize=(7.05, 1.90))
    fig.subplots_adjust(left=0.07, right=0.995, top=0.78, bottom=0.28)

    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x - width / 2, ncd_time, width=width, color="0.78", hatch="//",
           edgecolor="black", linewidth=0.5, label="CND")
    ax.bar(x + width / 2, reg_time, width=width, color="0.08",
           edgecolor="black", linewidth=0.5, label="RegND*")
    ax.set_yscale("log")
    ax.set_ylim(min(reg_time) * 0.55, max(ncd_time) * 7.5)
    ax.set_ylabel("runtime (s)", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.3)
    for i, sp in enumerate(speedups):
        ax.text(i, max(ncd_time[i], reg_time[i]) * 1.35, fmt_speedup(sp),
                ha="center", va="bottom", fontsize=7.0)
    ax.legend(loc="upper center", ncol=2, frameon=False, fontsize=7.5,
              bbox_to_anchor=(0.5, 1.16))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)

    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[runtime ] -> {out_path}")


def make_memory_compare(rows: list[dict], out_path: Path) -> None:
    """Memory comparison on the same hard matched cells."""
    data = collect_hard_cases(rows)
    labels = data["labels"]
    ncd_memory = data["ncd_memory"]
    reg_memory = data["reg_memory"]
    memory_ratios = data["memory_ratios"]

    fig, ax = plt.subplots(1, 1, figsize=(3.45, 2.25))
    fig.subplots_adjust(left=0.16, right=0.99, top=0.82, bottom=0.24)

    x = np.arange(len(labels))
    width = 0.36
    ax.bar(x - width / 2, ncd_memory, width=width, color="0.78", hatch="//",
           edgecolor="black", linewidth=0.5, label="CND")
    ax.bar(x + width / 2, reg_memory, width=width, color="0.08",
           edgecolor="black", linewidth=0.5, label="RegND*")
    ax.set_yscale("log")
    ax.set_ylim(min(reg_memory) * 0.55, max(ncd_memory) * 4.0)
    ax.set_ylabel("memory (GB)", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    for i, ratio in enumerate(memory_ratios):
        ax.text(i, max(ncd_memory[i], reg_memory[i]) * 1.30, fmt_speedup(ratio),
                ha="center", va="bottom", fontsize=6.5)
    ax.legend(loc="upper center", ncol=2, frameon=False, fontsize=7.5,
              bbox_to_anchor=(0.5, 1.15))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)

    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[memory  ] -> {out_path}")


def make_coverage(rows: list[dict], out_path: Path) -> None:
    """Completed-cell coverage under the same time and memory budgets."""
    coverage_labels = [GRAPH_SHORT[g] for g in DENSE_COVERAGE_GRAPHS]
    coverage = {"CND": [], "RegND*": []}
    coverage_text = {"CND": [], "RegND*": []}
    for graph in DENSE_COVERAGE_GRAPHS:
        for algo, label in [("REF", "CND"), ("RegNDC", "RegND*")]:
            rows_ga = [row for row in rows if row["graph"] == graph and row["algo"] == algo]
            ok = sum(1 for row in rows_ga if row["status"] == "OK")
            total = len(rows_ga)
            pct = 100.0 * ok / total if total else 0.0
            coverage[label].append(pct)
            coverage_text[label].append(f"{ok:,}/{total:,}")

    fig, ax = plt.subplots(1, 1, figsize=(3.45, 2.20))
    fig.subplots_adjust(left=0.17, right=0.99, top=0.82, bottom=0.22)

    width = 0.36
    cx = np.arange(len(DENSE_COVERAGE_GRAPHS))
    ax.bar(cx - width / 2, coverage["CND"], width=width, color="0.78",
           hatch="//", edgecolor="black", linewidth=0.5, label="CND")
    ax.bar(cx + width / 2, coverage["RegND*"], width=width, color="0.08",
           edgecolor="black", linewidth=0.5, label="RegND*")
    ax.set_ylim(0, 108)
    ax.set_ylabel("completed cells (%)", fontsize=8)
    ax.set_xticks(cx)
    ax.set_xticklabels(coverage_labels, fontsize=7)
    for i in range(len(DENSE_COVERAGE_GRAPHS)):
        ax.text(cx[i] - width / 2, max(coverage["CND"][i] + 3, 4),
                coverage_text["CND"][i], ha="center", va="bottom",
                rotation=90, fontsize=5.8)
        ax.text(cx[i] + width / 2, min(coverage["RegND*"][i] - 4, 96),
                coverage_text["RegND*"][i], ha="center", va="top",
                rotation=90, fontsize=5.8, color="white")
    ax.legend(loc="upper center", ncol=2, frameon=False, fontsize=7.5,
              bbox_to_anchor=(0.5, 1.15))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)

    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[coverage] -> {out_path}")


# ---------------------------------------------------------------------------
def make_regnd_scale(rows: list[dict], out_path: Path) -> None:
    """RegND-only runtime scaling with s and with r."""
    idx = index_rows(rows)
    dense_graphs = ["dblp-core30", "com-dblp", "web-it-2004"]
    styles = {
        "dblp-core30": dict(color="black", linestyle="-", marker="o"),
        "com-dblp": dict(color="0.35", linestyle="--", marker="s"),
        "web-it-2004": dict(color="0.05", linestyle=":", marker="^"),
    }

    fig, axes = plt.subplots(
        1, 2, figsize=(7.05, 2.05),
        gridspec_kw=dict(wspace=0.34),
    )
    fig.subplots_adjust(left=0.08, right=0.99, top=0.78, bottom=0.25)

    ax = axes[0]
    for graph in dense_graphs:
        pts = []
        for (gg, r, s, algo), row in idx.items():
            if gg == graph and algo == "RegNDC" and r == 3 and row["status"] == "OK":
                t = ok_time_seconds(row)
                if t is not None:
                    pts.append((s, t))
        pts.sort()
        if not pts:
            continue
        checkpoints = (4, 8, 16, 32, 64, 114, 224, 432)
        plot_pts = [(x, y) for x, y in pts if x in checkpoints]
        if len(plot_pts) < 2:
            step = max(1, len(pts) // 8)
            plot_pts = pts[::step]
        xs = [p[0] for p in plot_pts]
        ys = [p[1] for p in plot_pts]
        ax.plot(xs, ys, linewidth=1.35, markersize=3.4,
                markerfacecolor="white", markeredgewidth=0.8,
                label=GRAPH_SHORT[graph], **styles[graph])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$s$ at fixed $r=3$", fontsize=9)
    ax.set_ylabel("runtime (s)", fontsize=9)
    ax.set_title("(a) scaling with $s$", fontsize=9)
    ax.set_xticks([4, 8, 16, 32, 64, 128, 256, 432])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{int(x)}" if x in [4, 8, 16, 32, 64, 128, 256, 432] else ""))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)

    ax = axes[1]
    r_values = range(3, 8)
    for graph in dense_graphs:
        pts = []
        for r in r_values:
            row = regnd_row(idx, graph, r, 10)
            t = ok_time_seconds(row)
            if t is not None:
                pts.append((r, t))
        if not pts:
            continue
        ax.plot([p[0] for p in pts], [p[1] for p in pts],
                linewidth=1.3, markersize=3.6, markerfacecolor="white",
                markeredgewidth=0.8, label=GRAPH_SHORT[graph],
                **styles[graph])
    ax.set_yscale("log")
    ax.set_xticks(list(r_values))
    ax.set_xlabel(r"$r$ at fixed $s=10$", fontsize=9)
    ax.set_title("(b) scaling with $r$", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3,
               frameon=False, fontsize=8.2, bbox_to_anchor=(0.52, 1.00))
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[scale   ] -> {out_path}")


# ---------------------------------------------------------------------------
def make_heatmap(rows: list[dict], out_path: Path) -> None:
    """6 graphs x 2 algos heatmap. Cell color = log10(wall_ms / 1000).

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
        figsize=(2.15 * n_graphs, 1.45 * n_algos + 0.55),
        gridspec_kw=dict(wspace=0.18, hspace=0.28),
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
                              labelpad=20, ha="right", va="center")
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
    cbar_ax = fig.add_axes([0.93, 0.20, 0.012, 0.62])
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
               ncol=4, fontsize=8, frameon=False,
               bbox_to_anchor=(0.45, 0.00))

    fig.subplots_adjust(left=0.10, right=0.91, top=0.91, bottom=0.16)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[heatmap] -> {out_path}")


# ---------------------------------------------------------------------------
def make_memory_ratio(rows: list[dict], out_path: Path) -> None:
    """Memory ratio M_REF / M_V3LM vs s, one curve per r in R_VALUES,
    one panel per paper-6 graph that has any matched cell."""
    idx = index_rows(rows)

    # Build per-graph data: data[graph][r] = list of (s, ratio).
    # Use RegNDC label (newer bench snapshot of the RegND* algorithm);
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
    if n == 5:
        fig = plt.figure(figsize=(7.2, 3.95))
        grid = fig.add_gridspec(
            2, 6, left=0.08, right=0.99, top=0.78, bottom=0.13,
            wspace=0.95, hspace=0.62,
        )
        axes_list = [
            fig.add_subplot(grid[0, 0:2]),
            fig.add_subplot(grid[0, 2:4]),
            fig.add_subplot(grid[0, 4:6]),
            fig.add_subplot(grid[1, 1:3]),
            fig.add_subplot(grid[1, 3:5]),
        ]
        cols = 3
    else:
        cols = min(n, 3)
        rows_g = math.ceil(n / cols)
        fig, axes = plt.subplots(
            rows_g, cols,
            figsize=(2.35 * cols, 1.82 * rows_g + 0.65),
            squeeze=False,
            gridspec_kw=dict(wspace=0.38, hspace=0.55),
        )
        fig.subplots_adjust(left=0.08, right=0.99, top=0.78, bottom=0.13)
        axes_list = [axes[i // cols][i % cols] for i in range(rows_g * cols)]

    for idx_g, g in enumerate(plotted_graphs):
        ax = axes_list[idx_g]
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
        ax.yaxis.set_major_formatter(FuncFormatter(compact_log_tick))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.get_offset_text().set_visible(False)
        ax.set_title(PAPER6_PLAIN[g], fontsize=10)
        ax.set_xlabel("s", fontsize=9)
        if idx_g in (0, 3):
            ax.set_ylabel(r"$M_{\mathrm{NuclearCD}} / M_{\mathrm{RegND}}$", fontsize=9)
        ax.tick_params(axis="both", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Hide unused panels
    for ax in axes_list[len(plotted_graphs):]:
        ax.set_visible(False)

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
               bbox_to_anchor=(0.5, 0.98))

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
    print(f"  matched (REF OK & RegND* OK) cells: {len(matched)}")
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
    p.add_argument("--only", nargs="*",
                   choices=["main", "time-advantage", "memory-advantage",
                            "runtime", "compare-memory", "coverage",
                            "scale", "heatmap", "memory"],
                   default=["time-advantage", "runtime",
                            "memory-advantage", "coverage"])
    p.add_argument("--no-summary", action="store_true")
    args = p.parse_args()

    if not args.csv.exists():
        print(f"missing CSV: {args.csv}", file=sys.stderr); sys.exit(1)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.csv)
    print(f"[load] {len(rows)} rows from {args.csv}")

    if "main" in args.only:
        make_main_results(rows, args.out_dir / "exp_main_results.pdf")
    if "time-advantage" in args.only:
        make_time_advantage(rows, args.out_dir / "exp_time_advantage.pdf")
    if "memory-advantage" in args.only:
        make_memory_advantage(rows, args.out_dir / "exp_memory_advantage.pdf")
    if "runtime" in args.only:
        make_runtime_compare(rows, args.out_dir / "exp_runtime.pdf")
    if "compare-memory" in args.only:
        make_memory_compare(rows, args.out_dir / "exp_memory_compare.pdf")
    if "coverage" in args.only:
        make_coverage(rows, args.out_dir / "exp_coverage.pdf")
    if "scale" in args.only:
        make_regnd_scale(rows, args.out_dir / "exp_scale_regnd.pdf")
    if "heatmap" in args.only:
        make_heatmap(rows, args.out_dir / "exp_heatmap.pdf")
    if "memory" in args.only:
        make_memory_ratio(rows, args.out_dir / "exp_memory_r.pdf")
    if not args.no_summary:
        numerical_summary(rows)


if __name__ == "__main__":
    main()
