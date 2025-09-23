import os
from typing import List, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

from python.experiment import make_plots, make_plots_bazel
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
    'CBS': '////',
    'CBS-noHi': '',
    'ARC': 'xx',
    'Nuclear-YES': '...',
    'Nuclear-NO': 'oo',
    'data3-YES': '...',
    'data3-NO': 'oo',
    'data3': '///',
}


COLOUR_MAP = {
    'CBS': 'gray',
    'CBS-noHi': 'BLACK',
    'ARC': 'w',
    'Nuclear-YES': 'w',
    'Nuclear-NO': 'w',
    'data3-YES': 'w',
    'data3-NO': 'w',
    'data3': 'w',
}

def _style_axes(ax):
    """Apply clean academic style to an axes."""
    # Ticks
    ax.tick_params(axis='both', which='both', direction='out', length=4, width=0.8)
    ax.minorticks_on()
    # x-axis: keep exactly one major tick per label; disable all minor ticks
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
    # Legend: tidy defaults if present later
    return ax

# ----- Output directory for comparison plots -----
COMPARE_OUT = "/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/"
os.makedirs(COMPARE_OUT, exist_ok=True)

EPS_T = 1e-1   # minimum positive value for time (seconds) on log scale
EPS_M = 1   # minimum positive value for memory (MB) on log scale
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
                _lbl.set_rotation(30)
                _lbl.set_ha('right')
            _style_axes(ax1)
            ax1.set_yscale('log')
            ax1.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)
            # add legend for time plot
            _h, _l = ax1.get_legend_handles_labels()
            if _l:
                ax1.legend(loc='upper left', ncol=min(3, len(_l)), frameon=True)
            # 30 degree for x

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


if __name__ == "__main__":
    # --- Read standardized CSVs directly (monochrome plotting) ---
    # data1 = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataNew") # CBS
    # data2 = make_plots_bazel.arcData() # ARC
    # data3 = arcData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/nucleus_experiment.log") # Nucleus
    # data4 = make_plots.myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataHi") # CBS
    # Prepare/standardize three tables
    # df_a = _prep_data1(data1)
    # df_b = _prep_data2(data2)
    # df_c = _prep_data3(data3)
    # df_d = _prep_data1(data4)
    # df_d["source"] = "CBS-noHi"

    # df_a.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data1.csv", index=False)
    # df_b.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data2.csv", index=False)
    # df_c.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data3.csv", index=False)
    # df_d.to_csv("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/data4.csv", index=False)


    DATA_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data"
    df_a = pd.read_csv(os.path.join(DATA_DIR, "data1.csv"))  # CBS
    df_b = pd.read_csv(os.path.join(DATA_DIR, "data2.csv"))  # ARC
    df_c = pd.read_csv(os.path.join(DATA_DIR, "data3.csv"))  # Nuclear (YES/NO variants already encoded in 'source')
    df_d = pd.read_csv(os.path.join(DATA_DIR, "data4.csv"))  # CBS-noHi

    # If r is 1 or 2 in df_a, treat those runs as CBS-noHi for plotting
    if 'source' not in df_a.columns:
        df_a['source'] = 'CBS'
    if 'r' in df_a.columns:
        try:
            df_a['r'] = df_a['r'].astype(int)
        except Exception:
            pass
        df_a.loc[df_a['r'].isin([1, 2]), 'source'] = 'CBS-noHi'


    df_c = df_c[df_c["source"] == "Nuclear-NO"]
    # Ensure types and minimal sanity
    for d in (df_a, df_b, df_c, df_d):
    # for d in (df_b, df_c, df_d):
        if 'r' in d.columns: d['r'] = d['r'].astype(int)
        if 's' in d.columns: d['s'] = d['s'].astype(int)
        # If time is missing or abnormal, keep as NaN so downstream EPS handles it
        if 'exit_status' in d.columns:
            d.loc[d['exit_status'].notna() & (d['exit_status'] != 0), 'total_sec'] = np.nan

    # --- Build nested dicts (dataset -> r -> s -> value) for each CSV ---
    def _best_per_drs(d: pd.DataFrame) -> pd.DataFrame:
        # keep the smallest non-null total_sec per (dataset,r,s); if all NaN, keep first
        if d.empty:
            return d
        dd = d.copy()
        dd['__rank__'] = dd['total_sec'].fillna(np.inf)
        idx = dd.groupby(['dataset_name','r','s'])['__rank__'].idxmin()
        out = dd.loc[idx].drop(columns='__rank__').reset_index(drop=True)
        return out

    def _to_nested_dict(d: pd.DataFrame) -> Dict:
        time_map: Dict[str, Dict[int, Dict[int, float]]] = {}
        mem_map: Dict[str, Dict[int, Dict[int, float]]] = {}
        for _, row in d.iterrows():
            ds = str(row['dataset_name'])
            r  = int(row['r'])
            s  = int(row['s'])
            t  = row['total_sec']
            # memory (MB) with fallbacks
            if 'max_rss_mb' in row.index and pd.notna(row['max_rss_mb']):
                m = float(row['max_rss_mb'])
            elif 'max_rss_kb' in row.index and pd.notna(row['max_rss_kb']):
                try:
                    m = float(row['max_rss_kb']) / 1024.0
                except Exception:
                    m = None
            else:
                m = None
            # time
            time_map.setdefault(ds, {}).setdefault(r, {})[s] = None if pd.isna(t) else float(t)
            # memory
            mem_map.setdefault(ds, {}).setdefault(r, {})[s] = None if pd.isna(m) else float(m)
        return {'time': time_map, 'memory': mem_map}

    # Reduce each CSV to best-per-(dataset,r,s) then export dicts
    df_a_best = _best_per_drs(df_a)
    df_b_best = _best_per_drs(df_b)
    df_c_best = _best_per_drs(df_c)
    df_d_best = _best_per_drs(df_d)

    data1 = _to_nested_dict(df_a_best)  # CBS
    data2 = _to_nested_dict(df_b_best)  # ARC
    data3 = _to_nested_dict(df_c_best)  # Nuclear (YES/NO already merged by best pick)
    data4 = _to_nested_dict(df_d_best)  # CBS-noHi

    # Combine, then for each (dataset, s, r, source) keep the non-null shortest time
    combined = pd.concat([df_a, df_d, df_b, df_c], ignore_index=True, sort=False)
    combined = _dedup_min_time_per_source(combined)

    # Generate plots: for each dataset and each R, Y is time/memory, X is s
    paths = _plot_compare_by_r(combined, COMPARE_OUT)

    # Write a simple preview index
    idx = pd.DataFrame([{"path": p} for p in paths])
    idx.to_csv(os.path.join(COMPARE_OUT, "_plots_index.csv"), index=False)
    with open(os.path.join(COMPARE_OUT, "_preview.html"), "w", encoding="utf-8") as f:
        f.write("<html><body><h2>Comparison (CBS / CBS-noHi / ARC / Nuclear): Time & Memory (X=s, grouped by r)</h2>\n")
        for p in paths:
            f.write(f"<div><img src='{os.path.basename(p)}' width='720'></div>\n")
        f.write("</body></html>")

    # Debug: print high-level sizes of dicts
    def _dict_size(d):
        return sum(len(d['time'].get(ds, {})) for ds in d['time'])
    print(f"Dicts built | data1(CBS):{_dict_size(data1)} data2(ARC):{_dict_size(data2)} data3(Nuclear):{_dict_size(data3)} data4(CBS-noHi):{_dict_size(data4)}")

    print(f"Wrote {len(paths)} plots to {COMPARE_OUT}")
    print("Done.")