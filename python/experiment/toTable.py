import os
from typing import List, Dict, Any, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import make_plots
import make_plots_bazel
from python.experiment.make_plots_nuclear import arcData

# ===== Academic plotting style =====
ACADEMIC_FONTSIZE_BASE = 12
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

def _style_axes(ax):
    """Apply clean academic style to an axes."""
    # Ticks
    ax.tick_params(axis='both', which='both', direction='out', length=4, width=0.8)
    ax.minorticks_on()
    # Grid on y, light dashed
    ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.4)
    # Spines: keep left/bottom, fade top/right
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.8)
    # Legend: tidy defaults if present later
    return ax

# ----- Output directory for comparison plots -----
COMPARE_OUT = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/image_compare"
os.makedirs(COMPARE_OUT, exist_ok=True)

EPS_T = 1e-3   # minimum positive value for time (seconds) on log scale
EPS_M = 1e-3   # minimum positive value for memory (MB) on log scale
TIME_YMAX_SEC = 6 * 3600  # 6 hours upper bound for time axis

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
    prefer_order = ["CBS", "ARC", "Nuclear-YES", "Nuclear-NO", "data3-YES", "data3-NO", "data3"]
    present = list(series.keys())
    ordered = [k for k in prefer_order if k in present] + [k for k in sorted(present) if k not in prefer_order]

    m = max(1, len(ordered))
    group_width = 0.8
    bar_w = group_width / m
    offsets = (np.arange(m) - (m - 1) / 2.0) * bar_w

    x_idx = np.arange(len(x_vals))
    for i, lab in enumerate(ordered):
        ys = series[lab]
        ax.bar(x_idx + offsets[i], ys, width=bar_w, label=lab, edgecolor='black', linewidth=0.6)

    ax.set_xticks(x_idx)
    ax.set_xticklabels([str(v) for v in x_vals])
    # Keep bars fully visible at both ends
    ax.set_xlim(-0.5, len(x_vals) - 0.5)

    ax.set_xlabel("r")  # 注意：在 _plot_compare_by_r 中 X 实为 s，但沿用原脚本接口不改签名
    ax.set_ylabel(y_label)
    ax.set_title(title)
    # ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.5)
    leg = ax.legend(frameon=False, handlelength=1.6, handletextpad=0.4, borderaxespad=0.4)


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
            fig1 = plt.figure(figsize=(6, 2))
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
            _style_axes(ax1)
            ax1.set_yscale('log')
            ax1.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)
            plt.tight_layout()
            p1 = os.path.join(out_dir, f"{dataset_name}_r{r_val}_time_by_s.png")
            fig1.savefig(p1, dpi=300, bbox_inches='tight')
            plt.close(fig1)
            paths.append(p1)

            # --- Memory figure (Y=log) ---
            fig2 = plt.figure(figsize=(6,2))
            ax2 = fig2.gca()
            series_m: Dict[str, List[float]] = {}
            for src in sorted(dd["source"].unique()):
                mm = dd[dd["source"] == src].set_index("s")["max_rss_mb"]
                ys = []
                for s in s_vals:
                    v = mm.get(s, np.nan)
                    try:
                        fv = float(v)
                    except Exception:
                        fv = np.nan
                    ys.append(fv if (pd.notna(fv) and fv > 0) else EPS_M)
                series_m[src] = ys
            _grouped_bar(ax2, s_vals, series_m, f"{dataset_name}  r={r_val}  (Memory)", "Max RSS (MB)")
            _style_axes(ax2)
            ax2.set_yscale('log')
            ax2.set_ylim(bottom=EPS_M)
            plt.tight_layout()
            p2 = os.path.join(out_dir, f"{dataset_name}_r{r_val}_mem_by_s.png")
            fig2.savefig(p2, dpi=300, bbox_inches='tight')
            plt.close(fig2)
            paths.append(p2)

    return paths


def _plot_compare_by_s(df: pd.DataFrame, out_dir: str) -> List[str]:
    if df.empty:
        return []
    paths: List[str] = []

    for dataset_name in sorted(df["dataset_name"].unique()):
        sub = df[df["dataset_name"] == dataset_name]
        for s_val in sorted(sub["s"].unique()):
            dd = sub[sub["s"] == s_val]

            r_vals = sorted(dd["r"].unique().tolist())
            if not r_vals:
                continue

            # --- Time figure ---
            fig1 = plt.figure(figsize=(6,2))
            ax1 = fig1.gca()
            series_t: Dict[str, List[float]] = {}
            for src in sorted(dd["source"].unique()):
                mm = dd[dd["source"] == src].set_index("r")["total_sec"]
                ys = [ float(mm.get(r, np.nan)) if pd.notna(mm.get(r, np.nan)) else 0.0 for r in r_vals ]
                series_t[src] = ys
            _grouped_bar(ax1, r_vals, series_t, f"{dataset_name}  s={s_val}  (Time)", "Total time (s)")
            _style_axes(ax1)
            ax1.set_yscale('log')
            ax1.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)
            plt.tight_layout()
            p1 = os.path.join(out_dir, f"{dataset_name}_s{s_val}_time_by_r.png")
            fig1.savefig(p1, dpi=300, bbox_inches='tight')
            plt.close(fig1)
            paths.append(p1)

            # --- Memory figure ---
            fig2 = plt.figure(figsize=(6,2))
            ax2 = fig2.gca()
            series_m: Dict[str, List[float]] = {}
            for src in sorted(dd["source"].unique()):
                mm = dd[dd["source"] == src].set_index("r")["max_rss_mb"]
                ys = [ float(mm.get(r, np.nan)) if pd.notna(mm.get(r, np.nan)) else 0.0 for r in r_vals ]
                series_m[src] = ys
            _grouped_bar(ax2, r_vals, series_m, f"{dataset_name}  s={s_val}  (Memory)", "Max RSS (MB)")
            _style_axes(ax2)
            ax2.set_yscale('log')
            ax2.set_ylim(bottom=EPS_M)
            plt.tight_layout()
            p2 = os.path.join(out_dir, f"{dataset_name}_s{s_val}_mem_by_r.png")
            fig2.savefig(p2, dpi=300, bbox_inches='tight')
            plt.close(fig2)
            paths.append(p2)

    return paths


if __name__ == "__main__":
    # --- Use your existing loaders (DO NOT re-parse logs) ---
    data1 = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataNew") # CBS
    data2 = make_plots_bazel.arcData() # ARC
    data3 = arcData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/nucleus_experiment.log") # Nucleus
    data4 = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataHi") # CBS
    # Prepare/standardize three tables
    df_a = _prep_data1(data1)
    df_b = _prep_data2(data2)
    df_c = _prep_data3(data3)
    df_d = _prep_data1(data4)
    df_d["source"] = "CBS-noHi"


    df_a.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data1.csv", index=False)
    df_b.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data2.csv", index=False)
    df_c.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data3.csv", index=False)
    df_d.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data4.csv", index=False)