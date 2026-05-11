#!/usr/bin/env python3
"""
Generate every experimental figure (and headline tables) for vldbNuclearR1,
case-studies excluded.

Figures emitted (paper_plots/):
    fig_exp_endtoend.pdf     §7.2  3-graph speedup + RSS panel
    fig_exp_time.pdf         §7.2  10-graph wall-clock vs s
    fig_exp_mem.pdf          §7.3  10-graph peak RSS vs s
    fig_phase_breakdown.pdf  §7.4  load/build/peel stack (Pure vs REF)
    fig_stress_time.pdf      §7.7  synthetic |V|=1000 dense, time
    fig_stress_mem.pdf       §7.7  synthetic |V|=1000 dense, RSS
    fig_par_scaling.pdf      §6.6  parallel construction thread scaling
    fig_friendster.pdf       §7.8  com-friendster 1.8B-edge wall + RSS

Tables emitted (paper_plots/, .tex + .md):
    tab_bd_time / tab_bd_mem / tab_par / tab_local

Data fetch:
    Defaults to "smart": rsync only files older than `--max-stale` (sec) on
    server side or missing locally.  --fetch forces, --no-fetch disables.
    Server pattern matches the one used in the bench scripts:
        ssh wenqianz@tods2 via z5286111@cse.unsw.edu.au jumphost.

Usage:
    python3 make_paper_plots.py                 # smart fetch + all plots
    python3 make_paper_plots.py --fetch         # always rsync
    python3 make_paper_plots.py --no-fetch      # local only
    python3 make_paper_plots.py --only fig_par_scaling fig_exp_time
"""
from __future__ import annotations
import argparse, csv, math, os, statistics, subprocess, sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# ---------------------------------------------------------------------------
# Paths and server config

ROOT          = Path(__file__).resolve().parent
DATA_DIR      = ROOT / "paper_data"
OUT_DIR       = ROOT / "paper_plots"
OUT_DIR.mkdir(exist_ok=True)

SSH_USER      = "wenqianz"
SSH_HOST      = "tods2.cse.unsw.edu.au"
SSH_JUMP_USER = "z5286111"
SSH_JUMP_HOST = "cse.unsw.edu.au"
REMOTE_ROOT   = f"~/nclique/paper_data"

# CSV files we depend on. Map local path → (remote relative path, optional).
CSV_FILES: dict[Path, tuple[str, bool]] = {
    DATA_DIR / "01_main_benchmark_762.csv":              ("01_main_benchmark_762.csv",       True),
    DATA_DIR / "01_main_benchmark_v3.csv":               ("01_main_benchmark_v3.csv",        True),
    DATA_DIR / "02_breakdown_summary.csv":               ("02_breakdown_summary.csv",        True),
    DATA_DIR / "02_breakdown_summary_v3.csv":            ("02_breakdown_summary_v3.csv",     True),
    DATA_DIR / "03_breakdown_median.csv":                ("03_breakdown_median.csv",         True),
    DATA_DIR / "14_scalability_webgoogle.csv":           ("14_scalability_webgoogle.csv",    True),
    DATA_DIR / "15_stress_synthetic_dense.csv":          ("15_stress_synthetic_dense.csv",   True),
    DATA_DIR / "15_stress_synthetic_dense_v3.csv":       ("15_stress_synthetic_dense_v3.csv", True),
    DATA_DIR / "bench_par_sdct.csv":                     ("bench_par_sdct.csv",              True),
    DATA_DIR / "bench_local_v4.csv":                     ("bench_local_v4.csv",              True),
    DATA_DIR / "friendster_billion/bench_billion.csv":   ("friendster_billion/bench_billion.csv", True),
}

# ---------------------------------------------------------------------------
# Fetch helpers

def _ssh_proxy_arg() -> str:
    """The single -e argument that rsync needs to bounce through the jumphost."""
    return (
        f"ssh -o StrictHostKeyChecking=no "
        f"-o ProxyCommand='ssh -o StrictHostKeyChecking=no "
        f"-W {SSH_HOST}:22 {SSH_JUMP_USER}@{SSH_JUMP_HOST}'"
    )

def _need_fetch(local: Path, max_stale_sec: float) -> bool:
    if not local.exists(): return True
    age = (Path(local).stat().st_mtime)
    import time
    return (time.time() - age) > max_stale_sec

def fetch_data(force: bool, skip: bool, max_stale_sec: float) -> None:
    if skip:
        print("[fetch] --no-fetch: using local data only.")
        return
    targets = []
    for local, (remote_rel, _) in CSV_FILES.items():
        if force or _need_fetch(local, max_stale_sec):
            targets.append((local, remote_rel))
    if not targets:
        print(f"[fetch] all CSVs fresh (≤ {max_stale_sec/60:.0f} min); skipping rsync.")
        return
    print(f"[fetch] rsyncing {len(targets)} files from {SSH_USER}@{SSH_HOST}")
    for local, remote_rel in targets:
        local.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            "rsync", "-avz", "--partial",
            "-e", _ssh_proxy_arg(),
            f"{SSH_USER}@{SSH_HOST}:{REMOTE_ROOT}/{remote_rel}",
            str(local),
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"[fetch]   ok   {remote_rel}")
        except subprocess.CalledProcessError as e:
            print(f"[fetch]   skip {remote_rel}  ({e.stderr.strip().splitlines()[-1] if e.stderr else 'rsync failed'})")

# ---------------------------------------------------------------------------
# CSV loaders

def load_csv(path: Path) -> list[dict]:
    if not path.exists():
        print(f"[warn] missing {path}; figures depending on it will be skipped.")
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def load_main_benchmark() -> list[dict]:
    """Returns the main r=1 benchmark rows in a canonical schema:
        graph, r, s, algorithm, time_ms, memory_kB, status

    Prefers V3 (01_main_benchmark_v3.csv); falls back to legacy
    01_main_benchmark_762.csv. The CSV `algorithm` column stores legacy
    bench labels (`Pure` for our event-driven algorithm, `REF_R1` for the
    NuclearCD baseline). This loader normalises both to the paper-facing
    names `SPIN*` and `CND` so downstream callers see one schema.
    """
    v3_path = DATA_DIR / "01_main_benchmark_v3.csv"
    out: list[dict] = []
    if v3_path.exists():
        for r in load_csv(v3_path):
            if r.get("status") != "OK": continue
            algo = r.get("algorithm", "")
            # CSV stores legacy bench labels ("Pure" for our event-driven
            # algorithm, "REF_R1" for the baseline). Plotter uses the paper
            # names "SPINSTAR" and "CND" everywhere downstream, so map here.
            if algo == "Pure":   algo = "SPINSTAR"
            elif algo == "REF_R1": algo = "CND"
            # Use peel-only time as the algorithm metric — both V3 and
            # REF share the same SDCT_Fused build phase, so the speedup
            # claim is about peel.  Falls back to took_ms (build+peel)
            # then wall_ms when peel_ms is missing (older CSVs).
            time_ms = (r.get("peel_ms") or r.get("took_ms")
                       or r.get("wall_ms") or "")
            mem     = r.get("time_max_rss_kB") or r.get("memory_kB") or ""
            out.append({
                "graph":     r.get("graph", ""),
                "r":         r.get("r", "1"),
                "s":         r.get("s", ""),
                "algorithm": algo,
                "time_ms":   time_ms,
                "memory_kB": mem,
                "status":    "OK",
            })
        if out:
            print(f"[csv] using {v3_path.name} (Pure / V3 SOTA, {len(out)} OK rows)")
            return out
        print(f"[csv] {v3_path.name} present but no OK rows — falling back")
    legacy = DATA_DIR / "01_main_benchmark_762.csv"
    rows = load_csv(legacy)
    # Legacy CSV has algo values "Ours_ST" / "REF_R1"; normalise to SPIN*/CND.
    for r in rows:
        a = r.get("algorithm", "")
        if a == "Ours_ST":   r["algorithm"] = "SPINSTAR"
        elif a == "REF_R1":  r["algorithm"] = "CND"
    if rows:
        print(f"[csv] using legacy {legacy.name} ({len(rows)} rows)")
    return rows

def median_by(rows: list[dict], group_keys: tuple[str,...], value_key: str,
              filter_fn=None) -> dict:
    """Group rows by group_keys; return median of value_key per group."""
    g: dict = defaultdict(list)
    for r in rows:
        if filter_fn and not filter_fn(r): continue
        try:
            v = float(r[value_key])
        except (ValueError, KeyError, TypeError):
            continue
        g[tuple(r[k] for k in group_keys)].append(v)
    return {k: statistics.median(v) for k, v in g.items() if v}

def gmean(xs: Iterable[float]) -> float:
    xs = [x for x in xs if x and x > 0]
    if not xs: return float("nan")
    return math.exp(sum(math.log(x) for x in xs) / len(xs))

# ---------------------------------------------------------------------------
# Plot styling

SPINSTAR_COLOR = "#1f4e9c"  # solid blue
CND_COLOR  = "#c0392b"  # dashed red
GRID_KW    = dict(which="both", color="#cccccc", linestyle=":", linewidth=0.5)

plt.rcParams.update({
    "font.family":   "DejaVu Sans",
    "font.size":     9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.4,
    "axes.linewidth":  0.6,
    "savefig.bbox":    "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype":    42,   # editable text in PDF
    "ps.fonttype":     42,
})

# Each graph displayed name (LaTeX-style \texttt).
GRAPH_DISPLAY = {
    "com-amazon.ungraph":       "com-amazon",
    "com-dblp":                 "com-dblp",
    "com-orkut":                "com-orkut",
    "com-youtube":              "com-youtube",
    "soc-pokec-relationships":  "soc-pokec",
    "twitter_combined":         "twitter",
    "web-Google":               "web-Google",
    "web-Stanford":             "web-Stanford",
    "web-it-2004":              "web-it-2004",
    "wiki-Talk":                "wiki-Talk",
}

def save(fig, name: str) -> None:
    pdf = OUT_DIR / f"{name}.pdf"
    png = OUT_DIR / f"{name}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=180)
    plt.close(fig)
    print(f"[plot]  wrote {pdf.relative_to(ROOT)} (+ .png)")

# ---------------------------------------------------------------------------
# Figure 1: fig_exp_endtoend — 3-graph speedup + RSS ratio panel

def fig_exp_endtoend() -> None:
    """3-panel speedup + RSS ratio on shared log y-axis.

    Per skill rules: monochrome (black/gray), legend on a dedicated top row,
    no gridlines, hidden top/right spines.
    """
    from matplotlib.lines import Line2D
    from matplotlib.gridspec import GridSpec
    rows = load_main_benchmark()
    if not rows: return
    GRAPHS = ["com-youtube", "web-Stanford", "web-it-2004"]

    t_med = median_by(rows, ("graph","s","algorithm"), "time_ms",
                      filter_fn=lambda r: r["status"]=="OK")
    m_med = median_by(rows, ("graph","s","algorithm"), "memory_kB",
                      filter_fn=lambda r: r["status"]=="OK")

    # figsize > rendered: LaTeX scales 9.0 → 7.0 ⇒ text ~7pt rendered
    fig = plt.figure(figsize=(9.0, 3.0))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[0.22, 1.0],
                  hspace=0.6, wspace=0.12,
                  top=0.86, bottom=0.20, left=0.09, right=0.98)

    # Top legend row
    ax_leg = fig.add_subplot(gs[0, :])
    ax_leg.axis("off")
    handles = [
        Line2D([0], [0], color="black", lw=1.6, marker="o", ms=4,
               mec="black", mfc="black", label="Speedup (CND / SPIN*)"),
        Line2D([0], [0], color="#808080", lw=1.4, ls="--", marker="s", ms=4,
               mec="#808080", mfc="white", mew=1.2, label="RSS ratio (CND / SPIN*)"),
    ]
    ax_leg.legend(handles=handles, loc="center", ncol=2, frameon=False,
                  fontsize=10, handlelength=3.2, columnspacing=2.5)

    speedup_pool, mem_pool = [], []
    axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    for i, (ax, g) in enumerate(zip(axes, GRAPHS)):
        s_set = sorted({int(k[1]) for k in t_med if k[0]==g})
        sp_xs, sp_ys, mr_xs, mr_ys = [], [], [], []
        for s in s_set:
            t_ours = t_med.get((g, str(s), "SPINSTAR"))
            t_ref  = t_med.get((g, str(s), "CND"))
            m_ours = m_med.get((g, str(s), "SPINSTAR"))
            m_ref  = m_med.get((g, str(s), "CND"))
            if t_ours and t_ref and t_ours > 0:
                sp_xs.append(s); sp_ys.append(t_ref / t_ours)
                speedup_pool.append(t_ref / t_ours)
            if m_ours and m_ref and m_ours > 0:
                mr_xs.append(s); mr_ys.append(m_ref / m_ours)
                mem_pool.append(m_ref / m_ours)
        ax.plot(sp_xs, sp_ys, color="black", marker="o", ms=4, lw=1.5,
                mec="black", mfc="black")
        ax.plot(mr_xs, mr_ys, color="#808080", marker="s", ms=4, lw=1.3,
                ls="--", mec="#808080", mfc="white", mew=1.2)
        # Break-even line at ratio=1
        ax.axhline(1.0, color="#bbbbbb", ls=":", lw=0.8)
        ax.text(0.02, 0.04, r"SPIN* beats CND above $1$",
                transform=ax.transAxes, fontsize=8, color="#666",
                style="italic")
        ax.set_xlabel(r"$s$", fontsize=10)
        if i == 0:
            ax.set_ylabel("ratio (CND / SPIN*)", fontsize=10)
        ax.set_yscale("log"); ax.set_xscale("log")
        ax.set_title(GRAPH_DISPLAY.get(g, g), fontsize=10.5, color="#222")
        # Skill: no grid, hide top/right spines
        ax.grid(False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=9, colors="#444")

    if speedup_pool:
        sup = (f"{len(speedup_pool)} matched cells   "
               f"speedup gmean$\\,{{=}}\\,${gmean(speedup_pool):.2f}$\\times$ "
               f"max$\\,{{=}}\\,${max(speedup_pool):.2f}$\\times$   "
               f"RSS ratio gmean$\\,{{=}}\\,${gmean(mem_pool):.2f}$\\times$ "
               f"max$\\,{{=}}\\,${max(mem_pool):.2f}$\\times$")
        fig.text(0.5, 0.99, sup, ha="center", va="top",
                 fontsize=9.5, color="#222")
    save(fig, "fig_exp_endtoend")

# ---------------------------------------------------------------------------
# Figure 2 + 3: fig_exp_time / fig_exp_mem  — 10-graph grid, ours vs ref

def _grid_plot(rows: list[dict], value_key: str, ylabel: str, fname: str,
               spin_rows: dict = None, spin_timeouts: set = None) -> None:
    """7-graph grid with shared legend on a top row, monochrome (black/gray)
    palette, no gridlines, minimal chart junk."""
    if not rows: return
    from matplotlib.lines import Line2D
    from matplotlib.gridspec import GridSpec

    med = median_by(rows, ("graph","s","algorithm"), value_key,
                    filter_fn=lambda r: r["status"]=="OK")
    # Count matched (g,s) cells per graph (both algos OK)
    matched_cells_per_graph: dict = {}
    for (g, s, a) in med:
        if a != "SPINSTAR": continue
        if (g, s, "CND") in med:
            matched_cells_per_graph[g] = matched_cells_per_graph.get(g, 0) + 1
    MIN_CELLS = 3
    graphs = {g for g, n in matched_cells_per_graph.items() if n >= MIN_CELLS}
    # Stable order: small/medium graphs first, large web graphs last
    order = [g for g in [
        "com-amazon.ungraph","com-dblp","com-youtube","com-orkut",
        "twitter_combined","wiki-Talk","soc-pokec-relationships",
        "web-Google","web-Stanford","web-it-2004"] if g in graphs]
    n = len(order)
    # Pick a balanced grid: prefer (rows, cols) so rows*cols == n exactly,
    # else minimize empty cells. For our common counts (8 -> 4x2, 10 -> 5x2).
    cols = 5 if n >= 10 else (4 if n >= 8 else 3)
    rows_n = math.ceil(n / cols)

    # Layout: 1 thin row for legend + rows_n data rows.
    # figsize tuned so the figure renders at natural size in a two-column
    # paper (\linewidth ~7in for figure*); each panel ends up ~1.5" wide.
    # figure* in SIGMOD acmart sigconf is ~7.0" wide. Match figsize so
    # \includegraphics[width=\linewidth]{...} renders at 1:1 (no font shrink).
    # figsize larger than rendered (\linewidth in figure* ≈ 7"). LaTeX scales
    # everything down ⇒ text in figure renders slightly smaller than body
    # (matplotlib's default sans glyphs look heavier than the paper's serif
    # at nominally same pt; the down-scale evens it out).
    PAGE_W = 9.0
    fig = plt.figure(figsize=(PAGE_W, 0.6 + 1.5 * rows_n))
    gs = GridSpec(rows_n + 1, cols, figure=fig,
                  height_ratios=[0.20] + [1.0] * rows_n,
                  hspace=1.30, wspace=0.45,
                  top=0.92, bottom=0.13, left=0.08, right=0.98)

    # Legend axis (top, spans all columns)
    ax_leg = fig.add_subplot(gs[0, :])
    ax_leg.axis("off")
    handles = [
        Line2D([0], [0], color="black", lw=1.6, marker="o", ms=4,
               mec="black", mfc="black", label="SPIN* (ours)"),
        Line2D([0], [0], color="#808080", lw=1.4, ls="--", marker="s", ms=4,
               mec="#808080", mfc="white", mew=1.2, label="CND"),
        Line2D([0], [0], color="#444444", lw=1.4, ls=":", marker="^", ms=4.5,
               mec="#444444", mfc="white", mew=1.2, label="SPIN (ours)"),
    ]
    if spin_timeouts:
        handles.append(
            Line2D([0], [0], color="#444444", lw=0, marker="x", ms=6, mew=1.6,
                   label="SPIN $\\geq 1$h"))
    ax_leg.legend(handles=handles, loc="center", ncol=len(handles), frameon=False,
                  fontsize=11, handlelength=3.2, columnspacing=1.8)

    # Data axes
    axes = []
    for idx in range(n):
        ax = fig.add_subplot(gs[1 + idx // cols, idx % cols])
        axes.append(ax)
        g = order[idx]
        s_set = sorted({int(k[1]) for k in med if k[0]==g})
        ours = [med.get((g, str(s), "SPINSTAR")) for s in s_set]
        ref  = [med.get((g, str(s), "CND"))  for s in s_set]
        sx_o = [s for s, v in zip(s_set, ours) if v]
        vy_o = [v for v in ours if v]
        sx_r = [s for s, v in zip(s_set, ref) if v]
        vy_r = [v for v in ref  if v]
        ax.plot(sx_o, vy_o, color="black", marker="o", ms=3.5, lw=1.4,
                mec="black", mfc="black")
        ax.plot(sx_r, vy_r, color="#808080", marker="s", ms=3.5, lw=1.2,
                ls="--", mec="#808080", mfc="white", mew=1.2)
        # SPIN (LocalH) OK cells — if data exists for this graph
        if spin_rows:
            sx_sp = []
            vy_sp = []
            for s in s_set:
                v = spin_rows.get((g, str(s)))
                if v is not None:
                    sx_sp.append(s)
                    vy_sp.append(v)
            if sx_sp:
                ax.plot(sx_sp, vy_sp, color="#444444", marker="^", ms=4, lw=1.2,
                        ls=":", mec="#444444", mfc="white", mew=1.2)
        # SPIN TIMEOUT cells: × at the 1h cap (3.6e6 ms) — only on time plot
        if spin_timeouts:
            TIMEOUT_MS = 3.6e6
            sx_to = sorted(int(s) for (g_, s) in spin_timeouts if g_ == g)
            if sx_to:
                ax.plot(sx_to, [TIMEOUT_MS] * len(sx_to), color="#444444",
                        marker="x", ms=6, mew=1.6, lw=0)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$s$", fontsize=9)
        if idx % cols == 0:
            ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(GRAPH_DISPLAY.get(g, g), fontsize=9, color="#222", pad=2)
        # Skill: drop gridlines, drop top/right spines
        ax.grid(False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=8, colors="#444")
        # If only a few s values (e.g. com-orkut s=2..4), force integer ticks
        # instead of log-scientific (2x10^0). Threshold ≤ 5 distinct values.
        all_x = sorted(set(sx_o + sx_r))
        if 0 < len(all_x) <= 5:
            ax.set_xticks(all_x)
            ax.xaxis.set_major_formatter(mtick.ScalarFormatter())
            ax.xaxis.set_minor_locator(mtick.NullLocator())

    save(fig, fname)

def _load_spin_rows(value_key: str) -> tuple[dict, set]:
    """Read SPIN (LocalH) rows from bench_local_v4.csv.
    Returns (ok_values, timeouts) where:
      - ok_values: {(graph, s_str): value} for OK cells
      - timeouts: {(graph, s_str)} for cells that hit the 1h wall (only
        meaningful for the time plot; the mem plot ignores this set).
    """
    ok_values, timeouts = {}, set()
    csv_path = DATA_DIR / "bench_local_v4.csv"
    if not csv_path.exists(): return ok_values, timeouts
    for r in load_csv(csv_path):
        status = r.get("status")
        if status == "OK":
            try:
                if value_key == "time_ms":
                    raw = r.get("took_ms") or r.get("wall_ms")
                    v = float(raw)
                else:  # memory_kB
                    v = float(r["time_max_rss_kB"])
                ok_values[(r["graph"], r["s"])] = v
            except (KeyError, ValueError, TypeError):
                continue
        elif status == "TIMEOUT":
            timeouts.add((r["graph"], r["s"]))
    return ok_values, timeouts


def fig_exp_time() -> None:
    rows = load_main_benchmark()
    spin_rows, spin_timeouts = _load_spin_rows("time_ms")
    _grid_plot(rows, "time_ms", "wall time (ms)", "fig_exp_time",
               spin_rows, spin_timeouts)

def fig_exp_mem() -> None:
    rows = load_main_benchmark()
    spin_rows, _ = _load_spin_rows("memory_kB")
    _grid_plot(rows, "memory_kB", "peak RSS (kB)", "fig_exp_mem", spin_rows)

# ---------------------------------------------------------------------------
# Figure 4: fig_phase_breakdown — load/build/peel stack (Pure vs REF)

def fig_phase_breakdown() -> None:
    """V3 vs REF phase stack.  Prefers the V3 CSV (algorithm column =
    Pure / REF_R1, with 3 runs each — needs median); falls back to the
    legacy 02_breakdown_summary.csv (algo column = ours / ref, single
    pre-aggregated row per cell)."""
    v3_path = DATA_DIR / "02_breakdown_summary_v3.csv"
    legacy_path = DATA_DIR / "02_breakdown_summary.csv"
    cells = []
    if v3_path.exists():
        rows_v3 = load_csv(v3_path)
        # Aggregate to median over runs per (graph, s, algo).
        g: dict = defaultdict(lambda: defaultdict(list))
        for r in rows_v3:
            if r.get("status") != "OK": continue
            algo = "ours" if r.get("algorithm") == "Pure" else "ref"
            try:
                key = (r["graph"], int(r["s"]), algo)
                for col in ("load_ms", "build_ms", "peel_ms"):
                    if r.get(col):
                        g[key][col].append(float(r[col]))
            except (KeyError, ValueError):
                continue
        for (graph, s, algo), phases in g.items():
            if all(k in phases and phases[k] for k in ("load_ms","build_ms","peel_ms")):
                cells.append((graph, s, algo,
                              [statistics.median(phases["load_ms"]),
                               statistics.median(phases["build_ms"]),
                               statistics.median(phases["peel_ms"])]))
        if cells:
            print(f"[csv] using {v3_path.name} (Pure / V3 SOTA, {len(cells)} cells)")
    if not cells:
        rows = load_csv(legacy_path)
        if rows:
            print(f"[csv] using legacy {legacy_path.name} (ours = ST)")
        for r in rows:
            try:
                cells.append((r["graph"], int(r["s"]), r["algo"],
                              [float(r["load_ms"]), float(r["build_ms"]), float(r["peel_ms"])]))
            except (KeyError, ValueError):
                continue
    if not cells: return
    graphs = ["com-youtube","web-Stanford","web-it-2004"]
    fig, axes = plt.subplots(1, len(graphs), figsize=(13, 3.0), sharey=False)
    if len(graphs) == 1: axes = [axes]
    for ax, g in zip(axes, graphs):
        ours = sorted([(s, vals) for (gg, s, algo, vals) in cells
                       if gg==g and algo=="ours"], key=lambda x: x[0])
        ref  = sorted([(s, vals) for (gg, s, algo, vals) in cells
                       if gg==g and algo=="ref"],  key=lambda x: x[0])
        if not ours: continue
        s_vals = [s for s, _ in ours]
        idx = np.arange(len(s_vals))
        bar_w = 0.35
        # ours stack
        load_o = np.array([v[0] for _, v in ours])
        build_o= np.array([v[1] for _, v in ours])
        peel_o = np.array([v[2] for _, v in ours])
        ax.bar(idx - bar_w/2, load_o,                width=bar_w,
               color="#cccccc", label="load")
        ax.bar(idx - bar_w/2, build_o, bottom=load_o, width=bar_w,
               color=SPINSTAR_COLOR, label="build (SPIN*)")
        ax.bar(idx - bar_w/2, peel_o,  bottom=load_o+build_o, width=bar_w,
               color="#5dade2", label="peel (SPIN*)")
        # ref stack (next to it)
        if ref:
            ref_idx = []
            load_r, build_r, peel_r = [], [], []
            for s_o, _ in ours:
                match = next((v for s_r, v in ref if s_r == s_o), None)
                if match:
                    ref_idx.append(s_vals.index(s_o))
                    load_r.append(match[0]); build_r.append(match[1]); peel_r.append(match[2])
            load_r = np.array(load_r); build_r = np.array(build_r); peel_r = np.array(peel_r)
            ref_idx = np.array(ref_idx)
            ax.bar(ref_idx + bar_w/2, load_r,                  width=bar_w,
                   color="#cccccc", hatch="///", edgecolor="white")
            ax.bar(ref_idx + bar_w/2, build_r, bottom=load_r,   width=bar_w,
                   color=CND_COLOR, hatch="///", edgecolor="white",
                   label="build (CND)")
            ax.bar(ref_idx + bar_w/2, peel_r,  bottom=load_r+build_r, width=bar_w,
                   color="#f1948a", hatch="///", edgecolor="white",
                   label="peel (CND)")
        ax.set_xticks(idx); ax.set_xticklabels(s_vals)
        ax.set_xlabel(r"$s$"); ax.set_ylabel("time (ms)")
        ax.set_yscale("log")
        ax.set_title(GRAPH_DISPLAY.get(g, g))
        ax.grid(**GRID_KW, axis="y")
        if ax is axes[0]:
            ax.legend(loc="upper right", frameon=False, ncol=2, fontsize=7)
    fig.tight_layout()
    save(fig, "fig_phase_breakdown")

# ---------------------------------------------------------------------------
# Figure 5 + 6: fig_stress_time / fig_stress_mem — synthetic |V|=1000

def _stress_plot(value_key: str, ylabel: str, fname: str) -> None:
    """Synthetic |V|=1000 stress test, redesigned per skill rules.
    One panel per s value, monochrome (black SPIN* + gray CND), legend on top."""
    from matplotlib.lines import Line2D
    from matplotlib.gridspec import GridSpec
    v3_path = DATA_DIR / "15_stress_synthetic_dense_v3.csv"
    rows = []
    if v3_path.exists():
        for r in load_csv(v3_path):
            if r.get("algorithm") == "Pure":
                r = dict(r); r["algorithm"] = "SPINSTAR"
            rows.append(r)
        if rows:
            print(f"[csv] using {v3_path.name} (Pure / V3 SOTA, {len(rows)} rows)")
    if not rows:
        rows = load_csv(DATA_DIR / "15_stress_synthetic_dense.csv")
        if rows:
            print(f"[csv] using legacy 15_stress_synthetic_dense.csv")
    if not rows: return
    g = defaultdict(list)
    for r in rows:
        if r["status"] != "OK": continue
        try:
            g[(float(r["density"]), int(r["s"]), r["algorithm"])].append(float(r[value_key]))
        except (KeyError, ValueError):
            continue
    if not g: return
    med = {k: statistics.median(v) for k, v in g.items()}
    densities = sorted({k[0] for k in med})
    s_vals    = sorted({k[1] for k in med})
    n = len(s_vals)

    # figsize > rendered (single column \linewidth ≈ 3.4"): scale down ⇒
    # smaller text in PDF.
    rows_n = math.ceil(n / 2)
    fig = plt.figure(figsize=(4.5, 0.6 + 1.3 * rows_n))
    gs = GridSpec(rows_n + 1, 2, figure=fig,
                  height_ratios=[0.22] + [1.0] * rows_n,
                  hspace=1.0, wspace=0.50,
                  top=0.90, bottom=0.15, left=0.20, right=0.97)

    # Top legend
    ax_leg = fig.add_subplot(gs[0, :]); ax_leg.axis("off")
    handles = [
        Line2D([0], [0], color="black", lw=1.6, marker="o", ms=4,
               mec="black", mfc="black", label="SPINSTAR"),
        Line2D([0], [0], color="#808080", lw=1.4, ls="--", marker="s", ms=4,
               mec="#808080", mfc="white", mew=1.2, label="CND"),
    ]
    ax_leg.legend(handles=handles, loc="center", ncol=2, frameon=False,
                  fontsize=9, handlelength=2.5, columnspacing=2.0)

    for i, s in enumerate(s_vals):
        ax = fig.add_subplot(gs[1 + i // 2, i % 2])
        ours = [med.get((d, s, "SPINSTAR")) for d in densities]
        ref  = [med.get((d, s, "CND"))  for d in densities]
        d_o = [d for d, v in zip(densities, ours) if v]
        v_o = [v for v in ours if v]
        d_r = [d for d, v in zip(densities, ref) if v]
        v_r = [v for v in ref  if v]
        ax.plot(d_o, v_o, color="black", marker="o", ms=4, lw=1.5,
                mec="black", mfc="black")
        ax.plot(d_r, v_r, color="#808080", marker="s", ms=4, lw=1.3,
                ls="--", mec="#808080", mfc="white", mew=1.2)
        ax.set_yscale("log")
        ax.set_xlabel("density", fontsize=10)
        if i == 0:
            ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(f"$s={s}$", fontsize=10.5, color="#222")
        ax.grid(False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="both", labelsize=9, colors="#444")

    save(fig, fname)

def fig_stress_time() -> None:
    _stress_plot("time_ms", "wall time (ms)", "fig_stress_time")

def fig_stress_mem() -> None:
    _stress_plot("memory_kB", "peak RSS (kB)", "fig_stress_mem")

# ---------------------------------------------------------------------------
# Figure 7: fig_par_scaling — SDCT_Par4 thread scaling

def fig_par_scaling() -> None:
    rows = load_csv(DATA_DIR / "bench_par_sdct.csv")
    if not rows: return
    g: dict = defaultdict(list)
    for r in rows:
        if r["status"] != "OK" or not r.get("build_ms"): continue
        try:
            k = (r["graph"], int(r["s"]), int(r["threads"]))
            g[k].append(int(r["build_ms"]))
        except (KeyError, ValueError):
            continue
    med = {k: int(statistics.median(v)) for k, v in g.items()}
    THREADS = [1, 2, 4, 8, 16, 24, 32, 48, 64]
    pairs = sorted({(g_, s_) for (g_, s_, _) in med})
    if not pairs: return

    # Collapse to one line per (graph,s) — speedup vs T
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.0))
    ax_sp, ax_t = axes
    cmap = plt.cm.tab20(np.linspace(0, 1, len(pairs)))
    for ci, (g_, s_) in enumerate(pairs):
        t1 = med.get((g_, s_, 1))
        if not t1: continue
        xs = []; sp = []; tt = []
        for T in THREADS:
            v = med.get((g_, s_, T))
            if v:
                xs.append(T); sp.append(t1/v); tt.append(v)
        label = f"{GRAPH_DISPLAY.get(g_, g_)} s={s_}"
        ax_sp.plot(xs, sp, color=cmap[ci], marker="o", ms=3, label=label)
        ax_t.plot(xs, tt,  color=cmap[ci], marker="o", ms=3, label=label)
    # ideal line
    ax_sp.plot([1, 64], [1, 64], color="black", ls=":", lw=0.6, label="ideal")
    ax_sp.set_xscale("log"); ax_sp.set_yscale("log")
    ax_sp.set_xticks(THREADS); ax_sp.set_xticklabels(THREADS)
    ax_sp.xaxis.set_major_formatter(mtick.ScalarFormatter())
    ax_sp.set_xlabel("threads $T$"); ax_sp.set_ylabel("speedup vs $T{=}1$")
    ax_sp.set_title("Parallel construction strong scaling")
    ax_sp.grid(**GRID_KW)
    ax_t.set_xscale("log"); ax_t.set_yscale("log")
    ax_t.set_xticks(THREADS); ax_t.set_xticklabels(THREADS)
    ax_t.xaxis.set_major_formatter(mtick.ScalarFormatter())
    ax_t.set_xlabel("threads $T$"); ax_t.set_ylabel("build time (ms)")
    ax_t.set_title("Build time vs threads")
    ax_t.grid(**GRID_KW)
    if len(pairs) <= 24:
        ax_sp.legend(loc="upper left", frameon=False, ncol=2, fontsize=6,
                     bbox_to_anchor=(1.02, 1.0))
    fig.tight_layout()
    save(fig, "fig_par_scaling")

# ---------------------------------------------------------------------------
# Figure 8: fig_friendster — billion-edge wall + RSS vs s

def fig_friendster() -> None:
    rows = load_csv(DATA_DIR / "friendster_billion" / "bench_billion.csv")
    if not rows: return
    pts = []
    for r in rows:
        if r.get("status") != "OK": continue
        try:
            pts.append((int(r["s"]), float(r["wall_sec"]),
                       (float(r["mem_kB"])/1024/1024) if r.get("mem_kB") else None))
        except (KeyError, ValueError):
            continue
    if not pts: return
    pts.sort()
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.2))
    s_vals  = [p[0] for p in pts]
    wall    = [p[1] for p in pts]
    rss_gb  = [p[2] for p in pts]
    axes[0].plot(s_vals, wall, color=SPINSTAR_COLOR, marker="o", ms=4)
    axes[0].set_xlabel(r"$s$"); axes[0].set_ylabel("wall time (s)")
    axes[0].set_title("com-friendster, SPIN*, T=24")
    axes[0].grid(**GRID_KW)
    if any(r is not None for r in rss_gb):
        s_r = [s for s, r in zip(s_vals, rss_gb) if r is not None]
        v_r = [r for r in rss_gb if r is not None]
        axes[1].plot(s_r, v_r, color=CND_COLOR, marker="s", ms=4)
        axes[1].set_xlabel(r"$s$"); axes[1].set_ylabel("peak RSS (GB)")
        axes[1].set_title("com-friendster peak memory")
        axes[1].grid(**GRID_KW)
    else:
        axes[1].set_visible(False)
    fig.tight_layout()
    save(fig, "fig_friendster")

# ---------------------------------------------------------------------------
# Tables (LaTeX + Markdown)

def _emit_table(name: str, header: list[str], align: str, body: list[list[str]]) -> None:
    md_lines = ["| " + " | ".join(header) + " |",
                "|" + "|".join("---:" if a in "rR" else "---" for a in align) + "|"]
    for row in body:
        md_lines.append("| " + " | ".join(row) + " |")
    (OUT_DIR / f"{name}.md").write_text("\n".join(md_lines) + "\n")

    tex_lines = [r"\begin{tabular}{" + align + "}", r"\toprule",
                 " & ".join(header) + r" \\", r"\midrule"]
    for row in body:
        tex_lines.append(" & ".join(row) + r" \\")
    tex_lines += [r"\bottomrule", r"\end{tabular}"]
    (OUT_DIR / f"{name}.tex").write_text("\n".join(tex_lines) + "\n")
    print(f"[tab]   wrote {name}.tex (+ .md)")

def tab_par() -> None:
    rows = load_csv(DATA_DIR / "bench_par_sdct.csv")
    if not rows: return
    g: dict = defaultdict(list)
    for r in rows:
        if r["status"] != "OK" or not r.get("build_ms"): continue
        try:
            k = (r["graph"], int(r["s"]), int(r["threads"]))
            g[k].append(int(r["build_ms"]))
        except (KeyError, ValueError):
            continue
    med = {k: int(statistics.median(v)) for k, v in g.items()}
    THREADS = [1, 2, 4, 8, 16, 24, 32, 48, 64]
    pairs = sorted({(g_, s_) for (g_, s_, _) in med})
    body = []
    for g_, s_ in pairs:
        t1 = med.get((g_, s_, 1))
        if not t1: continue
        row = [GRAPH_DISPLAY.get(g_, g_), str(s_)]
        for T in THREADS:
            v = med.get((g_, s_, T))
            row.append(f"{v}" if v else "—")
        peak = max((t1/med[(g_, s_, T)] for T in THREADS[1:] if (g_, s_, T) in med),
                   default=None)
        row.append(f"{peak:.1f}×" if peak else "—")
        body.append(row)
    header = ["Graph", "$s$"] + [f"$T{{=}}{T}$" for T in THREADS] + ["best"]
    align  = "ll" + "r"*len(THREADS) + "r"
    _emit_table("tab_par", header, align, body)

def tab_local() -> None:
    local_rows = load_csv(DATA_DIR / "bench_local_v4.csv")
    main_rows = load_main_benchmark()
    if not local_rows or not main_rows: return

    pure = {
        (r["graph"], r["s"]): float(r["time_ms"])
        for r in main_rows
        if r.get("algorithm") == "SPINSTAR" and r.get("status") == "OK" and r.get("time_ms")
    }
    keep = {
        ("com-amazon.ungraph", "8"),
        ("com-dblp", "10"),
        ("com-youtube", "5"),
        ("soc-pokec-relationships", "5"),
        ("wiki-Talk", "5"),
        ("web-Google", "8"),
        ("web-Stanford", "10"),
    }
    body = []
    for r in local_rows:
        key = (r.get("graph"), r.get("s"))
        if key not in keep: continue
        graph, s = key
        status = r.get("status", "")
        p_ms = pure.get(key)
        if status == "OK":
            wall_s = float(r["wall_ms"]) / 1000.0
            rss_gb = float(r["time_max_rss_kB"]) / 1024.0 / 1024.0
            ratio = wall_s * 1000.0 / p_ms if p_ms else None
            body.append([
                GRAPH_DISPLAY.get(graph, graph), s, status,
                f"{wall_s:.1f}", f"{p_ms:.1f}" if p_ms else "—",
                f"{ratio:.0f}×" if ratio else "—", f"{rss_gb:.1f}",
            ])
        else:
            body.append([
                GRAPH_DISPLAY.get(graph, graph), s, status,
                "≥1800", f"{p_ms:.1f}" if p_ms else "—", "—", "—",
            ])
    order = {k: i for i, k in enumerate([
        "com-amazon", "com-dblp", "com-youtube", "soc-pokec",
        "wiki-Talk", "web-Google", "web-Stanford",
    ])}
    body.sort(key=lambda row: order.get(row[0], 999))
    header = ["Graph", "$s$", "SPIN", "SPIN wall (s)", "SPIN* wall (ms)", "ratio", "SPIN RSS (GB)"]
    align = "llrrrrr"
    _emit_table("tab_local", header, align, body)

def _load_breakdown_cells() -> dict:
    """Returns {(graph, s): {algo: (load, build, peel, mem_MB)}}, V3 first."""
    cells: dict = defaultdict(dict)
    v3_path = DATA_DIR / "02_breakdown_summary_v3.csv"
    if v3_path.exists():
        # Aggregate over runs: median per (graph, s, algo) per phase.
        agg: dict = defaultdict(lambda: defaultdict(list))
        for r in load_csv(v3_path):
            if r.get("status") != "OK": continue
            algo = "ours" if r.get("algorithm") == "Pure" else "ref"
            try:
                key = (r["graph"], int(r["s"]), algo)
                for col in ("load_ms", "build_ms", "peel_ms"):
                    if r.get(col):  agg[key][col].append(float(r[col]))
                if r.get("time_max_rss_kB"):
                    agg[key]["rss_mb"].append(float(r["time_max_rss_kB"]) / 1024.0)
                elif r.get("memory_kB"):
                    agg[key]["rss_mb"].append(float(r["memory_kB"]) / 1024.0)
            except (KeyError, ValueError):
                continue
        for (graph, s, algo), phases in agg.items():
            need = ("load_ms", "build_ms", "peel_ms")
            if all(k in phases and phases[k] for k in need):
                cells[(graph, s)][algo] = (
                    statistics.median(phases["load_ms"]),
                    statistics.median(phases["build_ms"]),
                    statistics.median(phases["peel_ms"]),
                    statistics.median(phases["rss_mb"]) if phases.get("rss_mb") else None,
                )
        if cells:
            print(f"[csv] tab_bd_*: using {v3_path.name} (Pure / V3 SOTA)")
            return cells
    # Legacy fallback (single pre-aggregated row, peak_rss_kb column).
    legacy = DATA_DIR / "02_breakdown_summary.csv"
    for r in load_csv(legacy):
        try:
            rss_mb = float(r["peak_rss_kb"]) / 1024.0 if r.get("peak_rss_kb") else None
            cells[(r["graph"], int(r["s"]))][r["algo"]] = (
                float(r["load_ms"]), float(r["build_ms"]), float(r["peel_ms"]), rss_mb)
        except (KeyError, ValueError):
            continue
    if cells:
        print(f"[csv] tab_bd_*: using legacy {legacy.name} (ours = ST)")
    return cells


def tab_bd_time() -> None:
    cells = _load_breakdown_cells()
    if not cells: return
    body = []
    for (g_, s_) in sorted(cells):
        ours = cells[(g_, s_)].get("ours")
        ref  = cells[(g_, s_)].get("ref")
        if not ours: continue
        row = [GRAPH_DISPLAY.get(g_, g_), str(s_)]
        # ours / ref are 4-tuples (load, build, peel, rss_mb); take the
        # first three for the time table.
        for v in ours[:3]: row.append(f"{v:.0f}")
        if ref:
            for v in ref[:3]: row.append(f"{v:.0f}")
        else:
            row += ["—", "—", "—"]
        body.append(row)
    header = ["Graph", "$s$",
              "SPIN* load", "SPIN* build", "SPIN* peel",
              "CND load", "CND build", "CND peel"]
    align  = "ll" + "r"*6
    _emit_table("tab_bd_time", header, align, body)


def tab_bd_mem() -> None:
    raw_cells = _load_breakdown_cells()
    if not raw_cells: return
    cells = defaultdict(dict)
    for (g_, s_), algos in raw_cells.items():
        for algo, vals in algos.items():
            if vals[3] is not None:
                cells[(g_, s_)][algo] = vals[3]
    body = []
    for (g_, s_) in sorted(cells):
        ours = cells[(g_, s_)].get("ours")
        ref  = cells[(g_, s_)].get("ref")
        if not (ours and ref): continue
        body.append([
            GRAPH_DISPLAY.get(g_, g_), str(s_),
            f"{ours:.0f}", f"{ref:.0f}", f"{ref/ours:.2f}×",
        ])
    header = ["Graph", "$s$", "SPIN* (MB)", "CND (MB)", "ratio"]
    align  = "llrrr"
    _emit_table("tab_bd_mem", header, align, body)

# ---------------------------------------------------------------------------
# Driver

PLOT_REGISTRY = {
    "fig_exp_endtoend":     fig_exp_endtoend,
    "fig_exp_time":         fig_exp_time,
    "fig_exp_mem":          fig_exp_mem,
    "fig_phase_breakdown":  fig_phase_breakdown,
    "fig_stress_time":      fig_stress_time,
    "fig_stress_mem":       fig_stress_mem,
    "fig_par_scaling":      fig_par_scaling,
    "fig_friendster":       fig_friendster,
    "tab_par":              tab_par,
    "tab_local":            tab_local,
    "tab_bd_time":          tab_bd_time,
    "tab_bd_mem":           tab_bd_mem,
}

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fetch",    action="store_true",
                   help="force rsync from server")
    p.add_argument("--no-fetch", action="store_true",
                   help="never contact server, use local CSVs only")
    p.add_argument("--max-stale", type=float, default=600,
                   help="seconds; if local CSV mtime is older than this, refetch (default 600)")
    p.add_argument("--only", nargs="+", choices=list(PLOT_REGISTRY),
                   help="produce only these named plots/tables")
    args = p.parse_args()

    fetch_data(force=args.fetch, skip=args.no_fetch, max_stale_sec=args.max_stale)

    targets = args.only if args.only else list(PLOT_REGISTRY)
    for name in targets:
        try:
            PLOT_REGISTRY[name]()
        except Exception as e:
            print(f"[err]  {name} failed: {e}")

    print(f"\nDone. Outputs in {OUT_DIR.relative_to(ROOT)}/")

if __name__ == "__main__":
    main()
