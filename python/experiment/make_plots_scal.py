import os
from typing import List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import make_plots, make_plots_bazel
from make_plots_nuclear import arcData

# ===== Academic plotting style =====
ACADEMIC_FONTSIZE_BASE = 18
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': ACADEMIC_FONTSIZE_BASE,
    'axes.labelsize': ACADEMIC_FONTSIZE_BASE,
    'axes.titlesize': ACADEMIC_FONTSIZE_BASE + 1,
    'legend.fontsize': ACADEMIC_FONTSIZE_BASE - 2,
    'xtick.labelsize': ACADEMIC_FONTSIZE_BASE - 2,
    'ytick.labelsize': ACADEMIC_FONTSIZE_BASE - 2,
    'savefig.dpi': 300,
    'figure.dpi': 120,
})

# Enforce monochrome look (no colors) and LaTeX‑friendly fonts
plt.rcParams.update({
    'mathtext.fontset': 'stix',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
plt.rcParams.update({
    'legend.edgecolor': 'black',
})
# Hatch patterns per series (stable order)
HATCH_MAP = {
    'CBS': '',
    'CBS-noHi': '///',
    'ARC': '',
    'Nuclear-YES': '',
    'Nuclear-NO': '///',
    'data3-YES': '',
    'data3-NO': '///',
    'data3': '...',
}

# Map to line style (instead of hatch)
LINESTYLE_MAP = {
    'CBS': 'solid', #'dashed','dotted'
    'CBS-noHi': 'solid',
    'ARC': 'solid',
    'Nuclear-YES': 'solid',
    'Nuclear-NO': 'solid',
    'data3-YES': 'solid',
    'data3-NO': 'solid',
    'data3': 'solid',
}

MARKER_MAP = {
    'CBS': 'o',
    'CBS-noHi': 'o',
    'ARC': 's',
    'Nuclear-YES': '^',
    'Nuclear-NO': '^',
    'data3-YES': 'X',
    'data3-NO': 'X',
    'data3': 'X',
}



COLOUR_MAP = {
    'CBS': 'black',  # 'lightgray','dimgray'
    'CBS-noHi': 'gray',
    'ARC': 'black',
    'Nuclear-YES': 'black',
    'Nuclear-NO': 'gray',
    'data3-YES': 'black',
    'data3-NO': 'gray',
    'data3': 'lightgray',
}

# ----- Legend helpers -----
_LEGEND_SAVED = False  # ensure we only write a standalone legend once per run

def _order_series_keys(present):
    prefer_order = ["CBS", "CBS-noHi", "ARC", "Nuclear-YES", "Nuclear-NO", "data3-YES", "data3-NO", "data3"]
    return [k for k in prefer_order if k in present] + [k for k in sorted(present) if k not in prefer_order]

def _save_standalone_legend(series_keys, out_dir, fname="legend_time_lines"):
    """
    Create a standalone legend figure (no axes) using line styles/markers consistent with _multi_line().
    Saves both PNG and EPS into out_dir.
    """
    if not series_keys:
        return
    ordered = _order_series_keys(list(series_keys))
    handles = []
    for lab in ordered:
        h = Line2D(
            [0], [0],
            color=COLOUR_MAP.get(lab, 'black'),
            linestyle=LINESTYLE_MAP.get(lab, 'solid'),
            marker=MARKER_MAP.get(lab, None),
            linewidth=1.6,
            markersize=4,
            label=lab,
        )
        handles.append(h)

    # Build a clean figure that only contains the legend
    fig = plt.figure(figsize=(6, 0.6 + 0.28 * ((len(ordered) - 1) // 4 + 1)))
    # Place legend at center; no axes
    leg = fig.legend(
        handles=handles,
        loc="center",
        ncol=min(4, len(ordered)),
        frameon=True,
        framealpha=1.0,     # avoid transparency warnings with PS/EPS backends
        edgecolor="black"
    )
    # Remove any axes by ensuring a blank canvas
    fig.gca().axis("off")

    png_path = os.path.join(out_dir, f"{fname}.png")
    # eps_path = os.path.join(out_dir, f"{fname}.eps")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    # fig.savefig(eps_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --- Standalone legend with fixed 5 algorithms, hatch-based ---
def _save_fixed_hatch_legend(out_dir, fname="legend_algorithms_fixed"):
    """
    Save a standalone legend using hatch-based bar patches for the fixed algorithms:
    CBS, CBS-noHi, ARC, Nuclear-YES, Nuclear-NO.
    """
    fixed = ["CBS", "CBS-noHi", "ARC", "Nuclear-YES", "Nuclear-NO"]
    handles = []
    for lab in fixed:
        handles.append(Patch(
            facecolor=COLOUR_MAP.get(lab, "white"),
            edgecolor="black",
            linewidth=0.8,
            hatch=HATCH_MAP.get(lab, ""),
            label=lab
        ))
    # Create a clean figure with only the legend
    fig = plt.figure(figsize=(6, 1.1))
    fig.legend(
        handles=handles,
        loc="center",
        ncol=min(5, len(handles)),
        frameon=True,
        framealpha=1.0,
        edgecolor="black"
    )
    fig.gca().axis("off")
    png_path = os.path.join(out_dir, f"{fname}.png")
    eps_path = os.path.join(out_dir, f"{fname}.eps")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(eps_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

# --- Standalone legend with fixed 5 algorithms, line/marker-based ---
def _save_fixed_line_legend(out_dir, fname="legend_algorithms_fixed_lines"):
    """
    Save a standalone legend using Line2D handles (line style + marker) for the fixed algorithms:
    CBS, CBS-noHi, ARC, Nuclear-YES, Nuclear-NO.
    """
    fixed = ["CBS", "CBS-noHi", "ARC", "Nuclear-YES", "Nuclear-NO"]
    handles = []
    for lab in fixed:
        handles.append(Line2D(
            [0], [0],
            color=COLOUR_MAP.get(lab, 'black'),
            linestyle=LINESTYLE_MAP.get(lab, 'solid'),
            marker=MARKER_MAP.get(lab, None),
            linewidth=1.6,
            markersize=4,
            label=lab
        ))
    fig = plt.figure(figsize=(6, 1.0))
    fig.legend(
        handles=handles,
        loc="center",
        ncol=min(5, len(handles)),
        frameon=True,
        framealpha=1.0,   # ensure no transparency for PS/EPS
        edgecolor="black"
    )
    fig.gca().axis("off")
    png_path = os.path.join(out_dir, f"{fname}.png")
    eps_path = os.path.join(out_dir, f"{fname}.eps")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(eps_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def _style_axes(ax):
    # """Apply clean academic style to an axes."""
    # # Ticks
    # ax.tick_params(axis='both', which='both', direction='out', length=4, width=0.8)
    # ax.minorticks_on()
    # # x-axis: keep exactly one major tick per label; disable all minor ticks
    # ax.tick_params(axis='x', which='major', length=3)
    # ax.tick_params(axis='x', which='minor', bottom=False)
    # ax.xaxis.set_minor_locator(NullLocator())
    # # Grid on y, light dashed
    # ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.4)
    # # Spines: keep left/bottom, fade top/right
    # for sp in ['top', 'right']:
    #     ax.spines[sp].set_visible(False)
    # for sp in ['left', 'bottom']:
    #     ax.spines[sp].set_linewidth(0.8)
    # # Legend: tidy defaults if present later
    # return ax
    # Ticks
    ax.tick_params(axis='both', which='both', direction='out', length=4, width=0.8)
    ax.minorticks_on()
    # x-axis: major ticks only, no minors
    ax.tick_params(axis='x', which='major', length=3)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.xaxis.set_minor_locator(NullLocator())
    # Grid on y, light dashed
    ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.4)
    # Spines: keep left/bottom, fade top/right
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.8)
    return ax

# ----- Output directory for comparison plots -----
COMPARE_OUT = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/"
os.makedirs(COMPARE_OUT, exist_ok=True)

EPS_T = 1e-1   # minimum positive value for time (seconds) on log scale
EPS_M = 1   # minimum positive value for memory (MB) on log scale
TIME_YMAX_SEC = 8 * 3600  # 6 hours upper bound for time axis

# --- Scalability by-R style (five lines = data sizes) ---
# --- Scalability by-R style (five lines = data sizes) ---
SIZE_ORDER = [20, 40, 60, 80, 100]
SIZE_LABEL = {20: "20%", 40: "40%", 60: "60%", 80: "80%", 100: "100%"}
SIZE_MARKER = {20: "o", 40: "s", 60: "^", 80: "X", 100: "D"}
# Distinct monochrome shades (light → dark; 100% darkest)
SIZE_COLOR = {
    20: "#666666",
    40: "#555555",
    60: "#444444",
    80: "#222222",
    100: "#000000",  # 最深
}
# Distinct dash styles per size (kept monochrome)
# Use either a named style or (offset, on_off_seq)
SIZE_LS = {
    20: "solid",           # 20%
    40: (0, (6, 3)),        # 40%
    60: (0, (4, 2, 1.5, 2)),# 60%  long-short pattern
    80: (0, (2, 2)),        # 80%  dotted-like
    100: (0, (8, 2, 2, 2)), # 100% long-short-short
}

# --- Only run scalability plotting (skip comparison/legend) ---
ONLY_SCALABILITY = True

# ===== Scalability helpers (merge dataScal + data4 to get sizes 20/40/60/80/100) =====
def _detect_size_and_base(df: pd.DataFrame, assume_100_if_no_suffix: bool = False) -> pd.DataFrame:
    d = df.copy()
    # Prefer explicit size columns if present
    explicit_size = None
    for c in ["size", "data_size", "scale", "edges_pct", "edges_percent", "percent", "Epercent", "ds", "ratio"]:
        if c in d.columns:
            try:
                explicit_size = pd.to_numeric(d[c].astype(str).str.extract(r"(\d+)")[0], errors="coerce")
            except Exception:
                explicit_size = pd.to_numeric(d[c], errors="coerce")
            break
    # Ensure required cols exist (include memory columns)
    for c in ["dataset_name", "r", "s", "total_sec", "tree_build_ms", "clique_index_ms", "prog_time_ms", "exit_status", "max_rss_mb", "max_rss_kb"]:
        if c not in d.columns:
            d[c] = np.nan
    # If total_sec missing, try to compute from *_ms
    if "total_sec" not in d.columns or d["total_sec"].isna().all():
        d["total_ms"] = d[["tree_build_ms", "clique_index_ms", "prog_time_ms"]].sum(axis=1, min_count=2)
        d["total_sec"] = d["total_ms"] / 1000.0
    # Abnormal exits -> NaN time
    d.loc[d["exit_status"].notna() & (d["exit_status"] != 0), "total_sec"] = np.nan
    # Memory MB: prefer max_rss_mb; else convert from max_rss_kb
    mem_mb = pd.to_numeric(d["max_rss_mb"], errors="coerce")
    if mem_mb.isna().all() and "max_rss_kb" in d.columns:
        mem_mb = pd.to_numeric(d["max_rss_kb"], errors="coerce") / 1024.0
    d["_mem_mb_"] = mem_mb

    # Parse size from dataset_name suffix: name-20 / name_20 / name(20) / name-20%
    name = d["dataset_name"].astype(str)
    m = name.str.extract(r"(?P<base>.*?)[\-\(_ ]?(?P<size>20|40|60|80|100)\s*%?p?\)?$")
    # remove everything after fitst .
    if explicit_size is not None:
        d["_size_"] = explicit_size
    else:
        d["_size_"] = pd.to_numeric(m["size"], errors="coerce")
    # Base dataset name without the size suffix
    base = m["base"].fillna(name)
    base = base.str.rstrip("-_ (")
    base = base.str.split(".").str[0]
    # Also drop trailing patterns such as ".p100", ".p80", ".p20", "-100", "_60" etc. from the base name
    base = base.str.replace(r"(?:[\._-]p?(?:20|40|60|80|100)\s*%?)$", "", regex=True)
    base = base.str.rstrip("._- ")
    d["_base_ds_"] = base
    # For the 100% dataset (data4), treat all rows as size=100 regardless of suffix
    if assume_100_if_no_suffix:
        d["_size_"] = 100
    # Drop rows still missing size
    d = d.dropna(subset=["_size_"])
    d["_size_"] = d["_size_"].astype(int)
    # Cast r,s to int
    d["r"] = pd.to_numeric(d["r"], errors="coerce").astype("Int64").fillna(0).astype(int)
    d["s"] = pd.to_numeric(d["s"], errors="coerce").astype("Int64").fillna(0).astype(int)
    return d[["_base_ds_", "_size_", "r", "s", "total_sec", "_mem_mb_"]].copy()

def  _build_scalability_dataframe(df_scal: pd.DataFrame, df4: pd.DataFrame) -> pd.DataFrame:
    """Return a cleaned dataframe with columns: _base_ds_, _size_ in {20,40,60,80,100}, r, s, total_sec."""
    a = _detect_size_and_base(df_scal, assume_100_if_no_suffix=False)
    b = _detect_size_and_base(df4, assume_100_if_no_suffix=True)
    d = pd.concat([a, b], ignore_index=True)
    d.to_csv(os.path.join(COMPARE_OUT, "scalability_raw_merged.csv"), index=False)
    d = d[d["_size_"].isin([20, 40, 60, 80, 100])]

    # For duplicates at same (_base_ds_, r, s, _size_), keep the smallest non-null time
    d["__rk__"] = d["total_sec"].fillna(np.inf)
    idx = d.groupby(["_base_ds_", "r", "s", "_size_"])["__rk__"].idxmin()
    d = d.loc[idx].drop(columns="__rk__").reset_index(drop=True)
    return d

def _plot_scalability_all(df: pd.DataFrame, out_dir: str) -> list:
    """
    For each dataset base and each (r,s), draw a line plot with X = data size (20/40/60/80/100),
    Y = total time (s, log). Save both PNG and EPS.
    Returns list of output paths.
    """
    if df.empty:
        return []
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ds in sorted(df["_base_ds_"].unique()):
        sub_ds = df[df["_base_ds_"] == ds]
        for (r_val, s_val), grp in sub_ds.groupby(["r", "s"]):
            xs = sorted(grp["_size_"].unique().tolist())
            # Build Y in the order of xs
            y = []
            for sz in xs:
                v = grp.loc[grp["_size_"] == sz, "total_sec"].astype(float)
                y.append(float(v.iloc[0]) if len(v) else np.nan)

            fig = plt.figure(figsize=(4.5, 3.5))
            ax = fig.gca()
            ax.plot(
                xs,
                y,
                label=f"{ds}  r={r_val} s={s_val}",
                color="black",
                linestyle="solid",
                marker="o",
                linewidth=1.6,
                markersize=4,
            )
            ax.set_xlabel("Data size (%)")
            ax.set_ylabel("Total time (s)")
            ax.set_xticks(xs)
            # rotate a little for dense ticks
            for _lbl in ax.get_xticklabels():
                _lbl.set_rotation(30)
                _lbl.set_ha('right')
            # ax.set_yscale("log")
            _style_axes(ax)
            ax.legend(loc="upper left", frameon=True, ncol=1)
            plt.tight_layout()

            base_fname = f"{ds.replace('/', '_').replace(' ', '_')}_r{int(r_val)}_s{int(s_val)}_scalability"
            png_path = os.path.join(out_dir, base_fname + ".png")
            eps_path = os.path.join(out_dir, base_fname + ".eps")
            fig.savefig(png_path, dpi=300, bbox_inches="tight")
            fig.savefig(eps_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            paths.extend([png_path, eps_path])
    return paths


# --- Scalability by-R, across s (five lines = sizes) ---
def _plot_scalability_by_r_across_s(df: pd.DataFrame, out_dir: str) -> list:
    """
    For each dataset base and each r, draw a figure with X = s (all s values),
    and five lines (sizes 20/40/60/80/100). Y = total time (s, log).
    """
    if df.empty:
        return []
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ds in sorted(df["_base_ds_"].unique()):
        sub_ds = df[df["_base_ds_"] == ds]
        for r_val in sorted(sub_ds["r"].unique()):
            grp = sub_ds[sub_ds["r"] == r_val]
            # collect all s across sizes
            s_vals = sorted(grp["s"].unique().tolist()) if "s" in grp.columns else []
            if not s_vals:
                # s not present in df -> cannot plot across-s lines
                continue
            fig = plt.figure(figsize=(4.5, 3.5))
            ax = fig.gca()
            for sz in SIZE_ORDER:
                gsz = grp[grp["_size_"] == sz]
                if gsz.empty:
                    continue
                mm = gsz.set_index("s")["total_sec"].astype(float)
                ys = []
                for s_ in s_vals:
                    v = mm.get(s_, np.nan)
                    ys.append(np.nan if pd.isna(v) else float(v))
                ax.plot(
                    s_vals,
                    ys,
                    label=SIZE_LABEL.get(sz, str(sz)),
                    color=SIZE_COLOR.get(sz, "black"),
                    linestyle=SIZE_LS.get(sz, "solid"),
                    marker=SIZE_MARKER.get(sz, "o"),
                    linewidth=1.8,
                    markersize=4.5,
                    markerfacecolor="none",
                    markeredgecolor=SIZE_COLOR.get(sz, "black"),
                )
            ax.set_xlabel("s")
            ax.set_ylabel("Total time (ss)")
            # rotate x tick labels for readability
            for _lbl in ax.get_xticklabels():
                _lbl.set_rotation(30)
                _lbl.set_ha('right')
            # ax.set_yscale("log")
            _style_axes(ax)
            # ax.legend(loc="upper left", frameon=True, framealpha=1.0)  # keep disabled unless needed
            plt.tight_layout()
            base_fname = f"{ds.replace('/', '_').replace(' ', '_')}_r{int(r_val)}_scalability_across_s"
            png_path = os.path.join(out_dir, base_fname + ".png")
            eps_path = os.path.join("/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/", base_fname + ".eps")

            fig.savefig(png_path, dpi=300, bbox_inches="tight")
            fig.savefig(eps_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            paths.extend([png_path, eps_path])
    return paths


# --- Memory scalability by-R, across s (five lines = sizes) ---
def _plot_mem_scalability_by_r_across_s(df: pd.DataFrame, out_dir: str) -> list:
    """
    For each dataset base and each r, draw a figure with X = s and five lines (20/40/60/80/100),
    Y = Max RSS in MB. Style matches time scalability plots.
    """
    if df.empty:
        return []
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for ds in sorted(df["_base_ds_"].unique()):
        sub_ds = df[df["_base_ds_"] == ds]
        for r_val in sorted(sub_ds["r"].unique()):
            grp = sub_ds[sub_ds["r"] == r_val]
            s_vals = sorted(grp["s"].unique().tolist()) if "s" in grp.columns else []
            if not s_vals:
                continue
            fig = plt.figure(figsize=(4.5, 3.5))
            ax = fig.gca()
            for sz in SIZE_ORDER:
                gsz = grp[grp["_size_"] == sz]
                if gsz.empty:
                    continue
                mm = gsz.set_index("s")["_mem_mb_"].astype(float) / 1024.0  # MB → GB
                ys = []
                for s_ in s_vals:
                    v = mm.get(s_, np.nan)
                    ys.append(np.nan if pd.isna(v) else float(v))
                ax.plot(
                    s_vals,
                    ys,
                    label=SIZE_LABEL.get(sz, str(sz)),
                    color=SIZE_COLOR.get(sz, "black"),
                    linestyle=SIZE_LS.get(sz, "solid"),
                    marker=SIZE_MARKER.get(sz, "o"),
                    linewidth=2.2,
                    markersize=5.0,
                    markerfacecolor="none",
                    markeredgecolor=SIZE_COLOR.get(sz, "black"),
                    markeredgewidth=1.2,
                )
            ax.set_xlabel("s")
            ax.set_ylabel("Memory (GB)")
            for _lbl in ax.get_xticklabels():
                _lbl.set_rotation(30)
                _lbl.set_ha('right')
            _style_axes(ax)
            # ax.legend(loc="upper left", frameon=True, framealpha=1.0)  # keep disabled unless needed
            plt.tight_layout()
            base_fname = f"{ds.replace('/', '_').replace(' ', '_')}_r{int(r_val)}_mem_scalability_across_s"
            png_path = os.path.join(out_dir, base_fname + ".png")
            eps_path = os.path.join("/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/", base_fname + ".eps")
            fig.savefig(png_path, dpi=300, bbox_inches="tight")
            fig.savefig(eps_path, dpi=300, bbox_inches="tight")
            plt.close(fig)
            paths.extend([png_path, eps_path])
    return paths

def _prep_data1(df: pd.DataFrame) -> pd.DataFrame:
    """Use existing DataFrame from make_plots.myData(); compute total_sec from three ms columns."""
    if df is None or df.empty:
        return pd.DataFrame(columns=["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"])

    d = df.copy()

    # Normalize column names we need
    need_cols = ["dataset_name", "s", "r", "tree_build_ms", "clique_index_ms", "prog_time_ms", "max_rss_kb", "exit_status"]
    for c in need_cols:
        if c not in d.columns:
            # create missing columns for robustness
            d[c] = np.nan

    # total time = sum of three ms fields; require at least two present (min_count=2) else NaN
    d["total_ms"] = d[["tree_build_ms", "clique_index_ms", "prog_time_ms"]].sum(axis=1, min_count=2)
    d["total_sec"] = d["total_ms"] / 1000.0

    # If abnormal exit, treat time as NaN
    d.loc[d["exit_status"].notna() & (d["exit_status"] != 0), "total_sec"] = np.nan

    # Memory: convert KB -> MB if available
    d["max_rss_mb"] = np.where(d.get("max_rss_mb") is not None and "max_rss_mb" in d.columns,
                               d.get("max_rss_mb"),
                               d["max_rss_kb"] / 1024.0)

    # Keep only columns we need
    d = d[["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status"]].copy()
    d["source"] = "CBS"

    # Drop rows missing essential keys
    d = d.dropna(subset=["dataset_name", "s", "r"])

    # Cast r,s to int where possible
    d["r"] = d["r"].astype(int)
    d["s"] = d["s"].astype(int)
    return d


def _prep_data2(df: pd.DataFrame) -> pd.DataFrame:
    """Use existing DataFrame from make_plots_bazel.arcData(); field total_sec already parsed."""
    if df is None or df.empty:
        return pd.DataFrame(columns=["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"])

    d = df.copy()
    for c in ["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status"]:
        if c not in d.columns:
            d[c] = np.nan

    # If abnormal exit, treat time as NaN
    d.loc[d["exit_status"].notna() & (d["exit_status"] != 0), "total_sec"] = np.nan

    # Keep only needed columns and add source
    d = d[["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status"]].copy()
    d["source"] = "ARC"

    d = d.dropna(subset=["dataset_name", "s", "r"])
    d["r"] = d["r"].astype(int)
    d["s"] = d["s"].astype(int)
    return d


def _prep_data3(df: pd.DataFrame) -> pd.DataFrame:
    """
    Use DataFrame from python.experiment.make_plots_nuclear.arcData();
    treat option (YES/NO) as different algorithms/series; keep only normal exits.
    - dataset_name: already cleaned in arcData
    - r, s: integers
    - total_sec: use wall_time_sec
    - memory: max_rss_mb
    - source: "data3-YES" or "data3-NO"
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"])

    d = df.copy()

    # Ensure columns exist for robustness
    for c in ["dataset_name", "s", "r", "user_time_sec", "max_rss_mb", "exit_status", "option"]:
        if c not in d.columns:
            d[c] = np.nan

    # Only keep normal exits for data3
    d = d[d["exit_status"].isin([0])].copy()

    # Time and memory
    d["total_sec"] = d["user_time_sec"]
    # Build source label by algorithm (option YES/NO -> different bars)
    d["source"] = d["option"].astype(str).map(lambda x: f"Nuclear-{x}")  # e.g., data3-YES, data3-NO

    # Keep only necessary columns
    d = d[["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"]].copy()

    # Drop rows missing essential keys and cast types
    d = d.dropna(subset=["dataset_name", "s", "r"])
    d["r"] = d["r"].astype(int)
    d["s"] = d["s"].astype(int)
    return d


def _dedup_min_time_per_source(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each (dataset_name, s, r, source), keep the row with the smallest non-null total_sec.
    If all rows in a group have NaN time, keep the first (needed for memory plots).
    Implemented via idxmin on a 'rank key' (NaN -> +inf) to avoid apply() warnings.
    """
    if df.empty:
        return df

    d = df.copy()
    d["__rank_time__"] = d["total_sec"].copy()
    d["__rank_time__"] = d["__rank_time__"].fillna(np.inf)

    keys = ["dataset_name", "s", "r", "source"]
    idx = d.groupby(keys)["__rank_time__"].idxmin()
    out = d.loc[idx].drop(columns="__rank_time__").sort_values(keys).reset_index(drop=True)
    return out


def _grouped_bar(ax, x_vals: List[int], series: Dict[str, List[float]], title: str, y_label: str):
    """
    Grouped bar helper with stable ordering and adaptive spacing.
    - Stable order: CBS, ARC, Nuclear-YES, Nuclear-NO, data3-YES, data3-NO, data3, then any others sorted.
    - Spacing: total group width fixed to 0.8; bar width = group_width / #series; centered offsets.
    """
    # Desired order to keep visual consistency across plots
    prefer_order = ["CBS","CBS-noHi", "ARC", "Nuclear-YES", "Nuclear-NO", "data3-YES", "data3-NO", "data3"]
    present = list(series.keys())
    ordered = [k for k in prefer_order if k in present] + [k for k in sorted(present) if k not in prefer_order]

    m = max(1, len(ordered))
    group_width = 0.8
    bar_w = group_width / m
    offsets = (np.arange(m) - (m - 1) / 2.0) * bar_w

    x_idx = np.arange(len(x_vals))
    for i, lab in enumerate(ordered):
        ys = series[lab]
        ax.bar(x_idx + offsets[i], ys, width=bar_w, label=lab, edgecolor='black', linewidth=0.6, facecolor=COLOUR_MAP.get(lab, ''), hatch=HATCH_MAP.get(lab, ''))

    ax.set_xticks(x_idx)
    ax.set_xticklabels([str(v) for v in x_vals])
    # Keep bars fully visible at both ends
    ax.set_xlim(-0.5, len(x_vals) - 0.5)

    ax.set_xlabel("r")  # 注意：在 _plot_compare_by_r 中 X 实为 s，但沿用原脚本接口不改签名
    ax.set_ylabel(y_label)
    # Title and legend removed for external control


def _multi_line(ax, x_vals: List[int], series: Dict[str, List[float]], title: str, y_label: str):
    """
    Multi-line helper with stable ordering.
    - Stable order: CBS, CBS-noHi, ARC, Nuclear-YES, Nuclear-NO, data3-YES, data3-NO, data3, then any others sorted.
    - No grouped offsets (lines share the same x positions).
    """
    present = list(series.keys())
    ordered = _order_series_keys(present)

    x_idx = np.arange(len(x_vals))

    for lab in ordered:
        ys = series[lab]
        ax.plot(
            x_idx,
            ys,
            label=lab,
            color=COLOUR_MAP.get(lab, 'black'),
            linestyle=LINESTYLE_MAP.get(lab, 'solid'),
            linewidth=1.6,
            marker=MARKER_MAP.get(lab, None),
            markersize=4,
        )

    ax.set_xticks(x_idx)
    ax.set_xticklabels([str(v) for v in x_vals])
    ax.set_xlim(-0.5, len(x_vals) - 0.5)

    ax.set_xlabel("r")  # (kept to match your original interface)
    ax.set_ylabel(y_label)


def _plot_compare_by_r(df: pd.DataFrame, out_dir: str) -> List[str]:
    """
    For each dataset and each r, draw two grouped-bar figures with X = s:
      (1) Time (Total time in seconds, log-scale)
      (2) Memory (Max RSS in MB, log-scale)
    Bars are grouped by 'source' (data1 vs data2).
    """
    if df.empty:
        return []
    paths: List[str] = []

    for dataset_name in sorted(df["dataset_name"].unique()):
        sub = df[df["dataset_name"] == dataset_name]
        for r_val in sorted(sub["r"].unique()):
            dd = sub[sub["r"] == r_val].copy()

            s_vals = sorted(dd["s"].unique().tolist())
            if not s_vals:
                continue

            # --- Time figure (Y=log) ---
            fig1 = plt.figure(figsize=(5,2))
            ax1 = fig1.gca()
            series_t: Dict[str, List[float]] = {}
            for src in sorted(dd["source"].unique()):
                mm = dd[dd["source"] == src].set_index("s")["total_sec"]
                ys = []
                for s in s_vals:
                    v = mm.get(s, np.nan)
                    try:
                        fv = float(v)
                    except Exception:
                        fv = np.nan
                    ys.append(fv if (pd.notna(fv) and fv > 0) else EPS_T)
                series_t[src] = ys
            _grouped_bar(ax1, s_vals, series_t, f"{dataset_name}  r={r_val}  (Time)", "Total time (s)")
            # rotate x tick labels by 30 degrees for readability
            for _lbl in ax1.get_xticklabels():
                _lbl.set_rotation(45)
                _lbl.set_ha('right')
            _style_axes(ax1)
            # ax1.set_yscale('log')
            ax1.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)
            # add legend for time plot
            _h, _l = ax1.get_legend_handles_labels()
            if _l:
                ax1.legend(loc='upper left', ncol=min(3, len(_l)), frameon=True)
            # 30 degree for x

            plt.tight_layout()
            p1 = os.path.join(out_dir, f"{dataset_name}_r{r_val}_time_by_s.png")
            fig1.savefig(p1, dpi=300, bbox_inches='tight')
            plt.close(fig1)
            paths.append(p1)


            # --- Time figure (Y=log, line chart) ---
            fig1 = plt.figure(figsize=(5, 2))
            ax1 = fig1.gca()
            series_t: Dict[str, List[float]] = {}
            na_pos: Dict[str, List[int]] = {}   # store indices of NaNs per series

            for src in sorted(dd["source"].unique()):
                mm = dd[dd["source"] == src].set_index("s")["total_sec"]

                ys = []
                na_idx = []
                for j, s in enumerate(s_vals):
                    v = mm.get(s, np.nan)
                    try:
                        fv = float(v)
                    except Exception:
                        fv = np.nan

                    if pd.isna(fv):
                        # NaN -> use TIME_YMAX_SEC and remember for red highlight
                        ys.append(TIME_YMAX_SEC)
                        na_idx.append(j)
                    else:
                        # keep previous behavior for non-positive values
                        ys.append(fv if fv > 0 else EPS_T)

                series_t[src] = ys
                na_pos[src] = na_idx

            # Save a standalone legend once (based on the currently present series)
            global _LEGEND_SAVED
            if not _LEGEND_SAVED:
                _save_standalone_legend(list(series_t.keys()), out_dir, fname="legend_time_lines")
                _LEGEND_SAVED = True

            _multi_line(ax1, s_vals, series_t, f"{dataset_name}  r={r_val}  (Time)", "Total time (s)")

            # overlay red markers where values were NaN
            x_idx = np.arange(len(s_vals))
            for lab, idxs in na_pos.items():
                if idxs:
                    ax1.scatter(
                        x_idx[idxs],
                        [TIME_YMAX_SEC] * len(idxs),
                        s=18,
                        color='red',
                        marker=MARKER_MAP.get(lab, 'o'),  # use same marker pattern
                        edgecolors='none',
                        zorder=4
                    )

            # rotate x tick labels for readability
            for _lbl in ax1.get_xticklabels():
                _lbl.set_rotation(45)
                _lbl.set_ha('right')

            _style_axes(ax1)
            # ax1.set_yscale('log')
            ax1.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)

            # add legend
            # _h, _l = ax1.get_legend_handles_labels()
            # if _l:
            #     ax1.legend(loc='upper left', ncol=min(3, len(_l)), frameon=True)

            plt.tight_layout()
            p1 = os.path.join(out_dir, f"{dataset_name}_r{r_val}_time_by_s.eps")
            fig1.savefig(p1, dpi=300, bbox_inches='tight')
            plt.close(fig1)
            paths.append(p1)




            # # --- Memory figure (Y=log) ---
            # fig2 = plt.figure(figsize=(6, 2))
            # ax2 = fig2.gca()
            # series_m: Dict[str, List[float]] = {}
            # for src in sorted(dd["source"].unique()):
            #     mm = dd[dd["source"] == src].set_index("s")["max_rss_mb"]
            #     ys = []
            #     for s in s_vals:
            #         v = mm.get(s, np.nan)
            #         try:
            #             fv = float(v)
            #         except Exception:
            #             fv = np.nan
            #         ys.append(fv if (pd.notna(fv) and fv > 0) else EPS_M)
            #     series_m[src] = ys
            # _grouped_bar(ax2, s_vals, series_m, f"{dataset_name}  r={r_val}  (Memory)", "Max RSS (MB)")
            # # rotate x tick labels by 30 degrees for readability
            # for _lbl in ax2.get_xticklabels():
            #     _lbl.set_rotation(30)
            #     _lbl.set_ha('right')
            # _style_axes(ax2)
            # ax2.set_yscale('log')
            # ax2.set_ylim(bottom=EPS_M)
            # # add legend for memory plot
            # _h, _l = ax2.get_legend_handles_labels()
            # if _l:
            #     ax2.legend(loc='upper left', ncol=min(3, len(_l)), frameon=True)
            # plt.tight_layout()
            # p2 = os.path.join(out_dir, f"{dataset_name}_r{r_val}_mem_by_s.png")
            # fig2.savefig(p2, dpi=300, bbox_inches='tight')
            # plt.close(fig2)
            # paths.append(p2)

    return paths


# if __name__ == "__main__":
    # # --- Read standardized CSVs directly (monochrome plotting) ---
    # # data1 = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataNew") # CBS
    # # data2 = make_plots_bazel.arcData() # ARC
    # # data3 = arcData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/nucleus_experiment.log") # Nucleus
    # # data4 = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataHi") # CBS
    #
    #
    # scal = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataScal") # CBS
    #
    # # df_a = _prep_data1(data1)
    # # df_b = _prep_data2(data2)
    # # df_c = _prep_data3(data3)
    # # df_d = _prep_data1(data4)
    # # df_d["source"] = "CBS-noHi"
    #
    #
    # df_scal = _prep_data1(scal)
    # df_scal["source"] = "Scalable-CBS"
    #
    # # df_a.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data1.csv", index=False)
    # # df_b.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data2.csv", index=False)
    # # df_c.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data3.csv", index=False)
    # # df_d.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data4.csv", index=False)
    # df_scal.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/dataScal.csv", index=False)
    #
    #
    # DATA_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/"
    # df_a = pd.read_csv(os.path.join(DATA_DIR, "data1.csv"))  # CBS
    # df_b = pd.read_csv(os.path.join(DATA_DIR, "data2.csv"))  # ARC
    # df_c = pd.read_csv(os.path.join(DATA_DIR, "data3.csv"))  # Nuclear (YES/NO variants already encoded in 'source')
    # df_d = pd.read_csv(os.path.join(DATA_DIR, "data4.csv"))  # CBS-noHi
    #
    # # If r is 1 or 2 in df_a, treat those runs as CBS-noHi for plotting
    # if 'source' not in df_a.columns:
    #     df_a['source'] = 'CBS'
    # if 'r' in df_a.columns:
    #     try:
    #         df_a['r'] = df_a['r'].astype(int)
    #     except Exception:
    #         pass
    #     df_a.loc[df_a['r'].isin([1, 2]), 'source'] = 'CBS-noHi'
    #
    #
    # df_c = df_c[df_c["source"] == "Nuclear-NO"]
    # # Ensure types and minimal sanity
    # for d in (df_a, df_b, df_c, df_d):
    #     # for d in (df_b, df_c, df_d):
    #     if 'r' in d.columns: d['r'] = d['r'].astype(int)
    #     if 's' in d.columns: d['s'] = d['s'].astype(int)
    #     # If time is missing or abnormal, keep as NaN so downstream EPS handles it
    #     if 'exit_status' in d.columns:
    #         d.loc[d['exit_status'].notna() & (d['exit_status'] != 0), 'total_sec'] = np.nan
    #
    # # --- Build nested dicts (dataset -> r -> s -> value) for each CSV ---
    # def _best_per_drs(d: pd.DataFrame) -> pd.DataFrame:
    #     # keep the smallest non-null total_sec per (dataset,r,s); if all NaN, keep first
    #     if d.empty:
    #         return d
    #     dd = d.copy()
    #     dd['__rank__'] = dd['total_sec'].fillna(np.inf)
    #     idx = dd.groupby(['dataset_name','r','s'])['__rank__'].idxmin()
    #     out = dd.loc[idx].drop(columns='__rank__').reset_index(drop=True)
    #     return out
    #
    # def _to_nested_dict(d: pd.DataFrame) -> Dict:
    #     time_map: Dict[str, Dict[int, Dict[int, float]]] = {}
    #     mem_map: Dict[str, Dict[int, Dict[int, float]]] = {}
    #     for _, row in d.iterrows():
    #         ds = str(row['dataset_name'])
    #         r  = int(row['r'])
    #         s  = int(row['s'])
    #         t  = row['total_sec']
    #         # memory (MB) with fallbacks
    #         if 'max_rss_mb' in row.index and pd.notna(row['max_rss_mb']):
    #             m = float(row['max_rss_mb'])
    #         elif 'max_rss_kb' in row.index and pd.notna(row['max_rss_kb']):
    #             try:
    #                 m = float(row['max_rss_kb']) / 1024.0
    #             except Exception:
    #                 m = None
    #         else:
    #             m = None
    #         # time
    #         time_map.setdefault(ds, {}).setdefault(r, {})[s] = None if pd.isna(t) else float(t)
    #         # memory
    #         mem_map.setdefault(ds, {}).setdefault(r, {})[s] = None if pd.isna(m) else float(m)
    #     return {'time': time_map, 'memory': mem_map}
    #
    # # Reduce each CSV to best-per-(dataset,r,s) then export dicts
    # df_a_best = _best_per_drs(df_a)
    # df_b_best = _best_per_drs(df_b)
    # df_c_best = _best_per_drs(df_c)
    # df_d_best = _best_per_drs(df_d)
    #
    # data1 = _to_nested_dict(df_a_best)  # CBS
    # data2 = _to_nested_dict(df_b_best)  # ARC
    # data3 = _to_nested_dict(df_c_best)  # Nuclear (YES/NO already merged by best pick)
    # data4 = _to_nested_dict(df_d_best)  # CBS-noHi
    #
    # # Combine, then for each (dataset, s, r, source) keep the non-null shortest time
    # combined = pd.concat([df_a, df_d, df_b, df_c], ignore_index=True, sort=False)
    # combined = _dedup_min_time_per_source(combined)
    #
    # if not ONLY_SCALABILITY:
    #     # --- Standalone legend (fixed 5 algorithms, hatch-based) ---
    #     _save_fixed_hatch_legend(COMPARE_OUT, fname="legend_algorithms_fixed")
    #     # --- Standalone legend (fixed 5 algorithms, line/marker-based) ---
    #     _save_fixed_line_legend(COMPARE_OUT, fname="legend_algorithms_fixed_lines")
    #
    #     # --- Standalone legend for ALL algorithms present across data ---
    #     try:
    #         all_series = sorted(set(combined['source'].dropna().astype(str).unique().tolist()))
    #     except Exception:
    #         all_series = []
    #     if all_series:
    #         _save_standalone_legend(all_series, COMPARE_OUT, fname="legend_all_algorithms")
    #         # prevent per-figure legend file generation
    #         _LEGEND_SAVED = True
    #
    #     # Generate plots: for each dataset and each R, Y is time/memory, X is s
    #     paths = _plot_compare_by_r(combined, COMPARE_OUT)
    #
    #     # Write a simple preview index
    #     idx = pd.DataFrame([{"path": p} for p in paths])
    #     idx.to_csv(os.path.join(COMPARE_OUT, "_plots_index.csv"), index=False)
    #     with open(os.path.join(COMPARE_OUT, "_preview.html"), "w", encoding="utf-8") as f:
    #         f.write("<html><body><h2>Comparison (CBS / CBS-noHi / ARC / Nuclear): Time & Memory (X=s, grouped by r)</h2>\n")
    #         for p in paths:
    #             f.write(f"<div><img src='{os.path.basename(p)}' width='720'></div>\n")
    #         f.write("</body></html>")
    #
    #     # Debug: print high-level sizes of dicts
    #     def _dict_size(d):
    #         return sum(len(d['time'].get(ds, {})) for ds in d['time'])
    #     print(f"Dicts built | data1(CBS):{_dict_size(data1)} data2(ARC):{_dict_size(data2)} data3(Nuclear):{_dict_size(data3)} data4(CBS-noHi):{_dict_size(data4)}")
    #
    #     print(f"Wrote {len(paths)} plots to {COMPARE_OUT}")
    #     print("Done.")
def _save_scalability_size_legend(out_dir, fname="legend_scal_sizes"):
    """Standalone legend for five data sizes (20/40/60/80/100),
    matching SIZE_COLOR / SIZE_LS / SIZE_MARKER. Writes PNG + EPS.
    """
    try:
        ordered = [sz for sz in SIZE_ORDER if sz in SIZE_LABEL]
        handles = []
        for sz in ordered:
            lab = SIZE_LABEL.get(sz, str(sz))
            h = Line2D(
                [0], [0],
                color=SIZE_COLOR.get(sz, 'black'),
                linestyle=SIZE_LS.get(sz, 'solid'),
                marker=SIZE_MARKER.get(sz, 'o'),
                linewidth=2.2,
                markersize=5.0,
                markerfacecolor='none',
                markeredgecolor=SIZE_COLOR.get(sz, 'black'),
                markeredgewidth=1.2,
                label=lab,
            )
            handles.append(h)
        if not handles:
            return
        fig = plt.figure(figsize=(4.8, 1.0))
        fig.legend(
            handles=handles,
            loc="center",
            ncol=len(handles),
            frameon=True,
            framealpha=1.0,   # 避免 EPS 透明度告警
            edgecolor="black",
        )
        fig.gca().axis('off')
        png_path = os.path.join(out_dir, f"{fname}.png")
        eps_path = os.path.join("/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/", f"{fname}.eps")
        fig.savefig(png_path, dpi=300, bbox_inches='tight')
        fig.savefig(eps_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
    except Exception as _e:
        print("[WARN] failed to save scalability size legend:", _e)
# ===== Scalability (ALL) =====
try:
    DATA_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/"
    df_scal_in = pd.read_csv(os.path.join(DATA_DIR, "dataScalNew.csv"))  # 20/40/60/80
    df_100_in = pd.read_csv(os.path.join(DATA_DIR, "dataScalNew.csv"))      # 100 (no suffix)
    # remove if dataset_name contains ".p"
    df_100_in = df_100_in[~df_100_in["dataset_name"].str.contains(r"\.p", na=False)].copy()
    # df_100_inb = pd.read_csv(os.path.join(DATA_DIR, "data1.csv"))      # 100 (no suffix)

    # Ensure 100% rows carry s values: left-join s grid from dataScal on (_base_ds_, r)
    df_scal_all = _build_scalability_dataframe(df_scal_in, df_100_in)
    # filter if s > 15
    df_scal_all = df_scal_all[df_scal_all["s"].le(15) | df_scal_all["s"].isna()].copy()
    df_scal_all = df_scal_all[df_scal_all["_base_ds_"].isin(["com-dblp","web-Google"])].copy()

    # df_scal_all.to_csv(os.path.join(COMPARE_OUT, "scalability_data_all.csv"), index=False)
    # tmp = pd.concat([df_scal_all, df_scal_in[["r", "s"]]], axis=1)
    s_grid = df_scal_all[["_base_ds_", "r", "s"]].dropna().drop_duplicates()
    mask_100 = (df_scal_all["_size_"] == 100) & (df_scal_all["s"].isna())
    # if mask_100.any():
    #     fill_100 = pd.merge(df_scal_all.loc[mask_100, ["_base_ds_", "r"]].drop_duplicates(), s_grid, on=["_base_ds_", "r"], how="left")
    #     # Expand rows for each missing (ds,r) to all s in grid and merge back times if uniquely keyed
    #     fill_100 = fill_100.dropna(subset=["s"]).drop_duplicates()
    #     df_scal_all = pd.concat([df_scal_all, fill_100.assign(_size_=100, total_sec=np.nan)], ignore_index=True).drop_duplicates(subset=["_base_ds_","r","s","_size_"], keep="first")

    SCAL_OUT = os.path.join(COMPARE_OUT, "scalability")
    os.makedirs(SCAL_OUT, exist_ok=True)
    scal_paths = _plot_scalability_by_r_across_s(df_scal_all, SCAL_OUT)
    pd.DataFrame({"path": scal_paths}).to_csv(os.path.join(SCAL_OUT, "_index_scalability_by_r_across_s.csv"), index=False)

    mem_paths = _plot_mem_scalability_by_r_across_s(df_scal_all, SCAL_OUT)
    pd.DataFrame({"path": mem_paths}).to_csv(os.path.join(SCAL_OUT, "_index_mem_scalability_by_r_across_s.csv"), index=False)

    _save_scalability_size_legend(SCAL_OUT, fname="legend_scal_sizes")
except Exception as e:
    print("[WARN] Scal®abi®lity plotting skipped due to error:", e)