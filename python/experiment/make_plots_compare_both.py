import os
import math
from typing import List, Dict, Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import NullLocator

import make_plots
import make_plots_bazel

# ----- Configuration -----
COMPARE_OUT = "/Users/zhangwenqian/Library/CloudStorage/Dropbox/应用/Overleaf/Nuclear CD/figure/"
DATA_DIR = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/"
REGENERATE_FROM_RAW = True  # Set True to re-run slow data extraction from raw logs

# --- Plotting Constants ---
EPS_T = 1e-1
EPS_M = 1
TIME_YMAX_SEC = 3600
MARKER_SIZE_PT = 4
TIMEOUT_MARKER_SIZE = int(MARKER_SIZE_PT * MARKER_SIZE_PT * 1.2)
YGRID_ON = False
HOLLOW_MARKERS = True
MARKER_EDGE_WIDTH = 0.9
MARKER_ALPHA = 1.0

# ===== Academic Plotting Style =====
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 12, 'axes.labelsize': 12,
    'axes.titlesize': 13, 'legend.fontsize': 10, 'xtick.labelsize': 10,
    'ytick.labelsize': 10, 'savefig.dpi': 300, 'figure.dpi': 120,
    'mathtext.fontset': 'stix', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'legend.edgecolor': 'black',
})

# --- Series Styling ---
HATCH_MAP = {
    'CBS': '', 'CBS-noHi': '///', 'ARB': '', 'ARB-16': '', 'Nuclear': '',
    'Nuclear-NO': '///', 'Nuclear-YES': '', 'data3-YES': '', 'data3-NO': '///', 'data3': '...',
    'ARB-noHi': '///', 'CBS3': '...',
}
LINESTYLE_MAP = {
    'CBS': 'solid', 'CBS-noHi': (0, (4, 1.5)), 'ARB': 'dashed', 'ARB-16': 'dashdot',
    'ARB-noHi': (0, (6, 2, 1.5, 2)), 'Nuclear': 'dotted', 'Nuclear-NO': (0, (2, 1)), 'Nuclear-YES': 'dotted',
    'data3-YES': (0, (8, 2)), 'data3-NO': (0, (3, 1, 1, 1)), 'data3': (0, (1, 1)),
    'CBS3': (0, (1, 1)),
}
MARKER_MAP = {
    'CBS': '', 'CBS-noHi': 'o', 'ARB': '', 'ARB-16': 'X', 'Nuclear': '',
    'Nuclear-NO': '^', 'Nuclear-YES': '', 'data3-YES': 'X', 'data3-NO': 'X', 'data3': 'X',
    'ARB-noHi': 's', 'CBS3': 'D',
}
COLOUR_MAP = {
    'CBS': 'black', 'CBS-noHi': 'lightgray', 'ARB': 'olive', 'ARB-16': 'black',
    'Nuclear': 'skyblue', 'Nuclear-NO': 'gray', 'Nuclear-YES': 'skyblue', 'data3-YES': 'black',
    'data3-NO': 'gray', 'data3': 'lightgray', 'ARB-noHi': 'gray',
    'CBS3': 'darkorange',
}


# --- Helper Functions ---

def _order_series_keys(present):
    prefer_order = ["CBS", "CBS-noHi", "CBS3", "ARB", "ARB-16", "ARB-noHi", "Nuclear", "Nuclear-YES", "Nuclear-NO", "data3-YES", "data3-NO",
                    "data3"]
    return [k for k in prefer_order if k in present] + [k for k in sorted(present) if k not in prefer_order]


def _get_smart_ticks(labels: List[Any], max_ticks: int = 10) -> (List[int], List[Any]):
    """
    Selects a reasonable number of ticks from a list of labels, ensuring the first
    and last are always included.
    """
    n = len(labels)
    if n <= max_ticks:
        return np.arange(n), labels

    # Generate evenly spaced indices, including the first and the last
    indices = np.linspace(0, n - 1, num=max_ticks)
    indices = np.unique(np.round(indices).astype(int))
    
    positions = indices
    display_labels = [labels[i] for i in indices]
    
    return positions, display_labels


def _prepare_cbs_data(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare CBS data from make_plots.myData()."""
    if df is None or df.empty: return pd.DataFrame()
    d = df.copy()
    for c in ["tree_build_ms", "clique_index_ms", "prog_time_ms", "max_rss_kb", "exit_status"]:
        if c not in d.columns: d[c] = np.nan
    d["total_sec"] = d[["tree_build_ms", "clique_index_ms", "prog_time_ms"]].sum(axis=1, min_count=2) / 1000.0
    d.loc[d["exit_status"] != 0, "total_sec"] = np.nan
    d["max_rss_mb"] = d["max_rss_kb"] / 1024.0
    d["source"] = "CBS"
    if 'r' in d.columns:
        # For r=1 and r=2, CBS is equivalent to CBS-noHi
        d.loc[d['r'].isin([1, 2]), 'source'] = 'CBS-noHi'
    return d[["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"]].dropna(
        subset=["dataset_name", "s", "r"])


def _prepare_arb_data(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare ARB data from make_plots_bazel.arcData()."""
    if df is None or df.empty: return pd.DataFrame()
    d = df.copy()
    d.loc[d["exit_status"] != 0, "total_sec"] = np.nan
    d["threads"] = pd.to_numeric(d.get("threads"), errors="coerce").fillna(1).astype(int)
    d["source"] = d["threads"].map(lambda x: "ARB" if x == 1 else f"ARB-{x}")
    return d[["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"]].dropna(
        subset=["dataset_name", "s", "r"])


def _prepare_nuclear_data(df: pd.DataFrame) -> pd.DataFrame:
    """Prepare Nuclear data from its specific log format."""
    if df is None or df.empty: return pd.DataFrame()
    d = df.copy()
    d = d[d["exit_status"] == 0]  # Only keep normal exits
    d["total_sec"] = d["user_time_sec"]
    d["source"] = "Nuclear-" + d["option"].astype(str)
    return d[["dataset_name", "s", "r", "total_sec", "max_rss_mb", "exit_status", "source"]].dropna(
        subset=["dataset_name", "s", "r"])


def _dedup_min_time_per_source(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates data to find the best (minimum) total_sec and max_rss_mb independently
    for each (dataset, s, r, source) group, considering only successful runs.
    """
    if df.empty: return df

    # Filter to ensure we only consider metrics from successful runs (exit_status == 0)
    # We use a copy to avoid modifying the original dataframe slice
    d = df.copy()

    # Set metrics of failed runs to NaN so they are ignored by min()
    # total_sec is likely already handled, but we do it for max_rss_mb too
    failed_mask = d['exit_status'] != 0
    d.loc[failed_mask, 'total_sec'] = np.nan
    d.loc[failed_mask, 'max_rss_mb'] = np.nan

    # Group and take the minimum (best) value for each metric independently
    # min() ignores NaNs. If all are NaN (all failed), result is NaN.
    return d.groupby(['dataset_name', 's', 'r', 'source'], as_index=False)[
        ['total_sec', 'max_rss_mb', 'exit_status']].min()


def _style_axes(ax):
    """Apply consistent styling to plot axes."""
    ax.tick_params(axis='both', which='major', direction='out', length=4, width=0.8)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())
    ax.grid(YGRID_ON)
    for spine in ['top', 'right']: ax.spines[spine].set_visible(False)
    for spine in ['left', 'bottom']: ax.spines[spine].set_linewidth(0.8)


def _save_standalone_legend(ax, out_dir: str, fname: str, ncol: int = None):
    """Extracts legend from axes and saves it as a standalone figure."""
    handles, labels = ax.get_legend_handles_labels()
    if not labels: return
    fig = plt.figure(figsize=(max(2.5, 0.7 * len(labels)), 0.9))
    fig.legend(handles, [lab.replace("CBS", "CND") for lab in labels], loc="center",
               ncol=(ncol or min(4, len(labels))), frameon=True, framealpha=1.0, edgecolor="black")
    fig.gca().axis("off")
    fig.savefig(os.path.join(out_dir, f"{fname}.eps"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# --- Analysis Function ---

def _analyze_cbs_vs_nohi_speedup(df: pd.DataFrame):
    """Analyzes and prints speedup comparison between CBS and CBS-noHi for each dataset."""
    if df.empty: return

    analysis_df = df[df["source"].isin(["CBS", "CBS-noHi"]) & (df["r"] > 2) & (df["s"] <= 15)].copy()
    analysis_df.loc[analysis_df["exit_status"] != 0, "total_sec"] = np.nan
    analysis_df = _dedup_min_time_per_source(analysis_df)

    for ds in sorted(analysis_df["dataset_name"].unique()):
        ds_df = analysis_df[analysis_df["dataset_name"] == ds].dropna(subset=['total_sec'])
        if ds_df.empty: continue

        piv = ds_df.pivot_table(index=["r", "s"], columns="source", values="total_sec", aggfunc="min").dropna()
        print(f"\n[Speedup Analysis] {ds}: Comparing {len(piv)} (r,s) pairs for CBS vs CBS-noHi")
        if piv.empty: continue

        for (r, s), row in piv.iterrows():
            if "CBS" in row and "CBS-noHi" in row and pd.notna(row['CBS']) and pd.notna(row['CBS-noHi']):
                pct = (row["CBS-noHi"] / row["CBS"] - 1.0) * 100.0
                print(f"  r={r}, s={s}: {pct:+.1f}% (CBS={row['CBS']:.3g}s, noHi={row['CBS-noHi']:.3g}s)")


# --- Main Plotting Functions ---

def _plot_cbs_vs_nohi_long_comparison(df: pd.DataFrame, out_dir: str):
    """Generates a 'long' bar plot with (r,s) pairs on the X-axis for CBS vs CBS-noHi."""
    if df.empty: return

    plot_df = df[df["source"].isin(["CBS", "CBS-noHi"]) & (df["r"] > 2) & (df["s"] <= 15)].copy()
    plot_df.loc[plot_df["exit_status"] != 0, "total_sec"] = np.nan
    plot_df = _dedup_min_time_per_source(plot_df)

    for ds in sorted(plot_df["dataset_name"].unique()):
        ds_df = plot_df[plot_df["dataset_name"] == ds]
        valid_df = ds_df.dropna(subset=['total_sec'])
        if valid_df.empty: continue

        pairs = sorted(valid_df[['r', 's']].drop_duplicates().itertuples(index=False, name=None))
        if not pairs: continue

        pivot = valid_df.pivot_table(index=['r', 's'], columns='source', values='total_sec', aggfunc='min')
        pivot = pivot.reindex(pd.MultiIndex.from_tuples(pairs, names=['r', 's']))
        series = {
            src: pivot[src].tolist()
            for src in ["CBS", "CBS-noHi"] if src in pivot.columns and not pivot[src].isnull().all()
        }
        if not series: continue

        fig, ax = plt.subplots(figsize=(5, 2))
        n_series = len(series)
        x_idx = np.arange(len(pairs))
        width = 0.85 / n_series
        offsets = (np.arange(n_series) - (n_series - 1) / 2.0) * width
        ordered_keys = [k for k in ["CBS", "CBS-noHi"] if k in series]
        for i, lab in enumerate(ordered_keys):
            ax.bar(x_idx + offsets[i], series[lab], width=width, label=lab, edgecolor='black', linewidth=0.6,
                   facecolor=COLOUR_MAP.get(lab, 'lightgray'), hatch=HATCH_MAP.get(lab, ''))

        ax.set_xticks(x_idx)
        ax.set_xticklabels([f"({r},{s})" for r, s in pairs], rotation=60, ha='right')
        ax.set_xlim(-0.5, len(pairs) - 0.5)
        ax.set_yscale("log")
        ax.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)
        ax.set_xlabel("(r, s)")
        ax.set_ylabel("Time (s)")
        _style_axes(ax)

        _save_standalone_legend(ax, out_dir, f"{ds}_legend_long", ncol=len(series))
        fig.savefig(os.path.join(out_dir, f"abc_comp_{ds}_time_by_rs_long.eps"), dpi=300, bbox_inches='tight')
        plt.close(fig)


def _plot_compare_by_r(df: pd.DataFrame, out_dir: str):
    """For each dataset and r, draw a line chart for time vs s across algorithms."""
    if df.empty: return

    for ds in sorted(df["dataset_name"].unique()):
        for r_val in sorted(df[df["dataset_name"] == ds]["r"].unique()):
            plot_df = df[(df["dataset_name"] == ds) & (df["r"] == r_val)]

            # Filter for r=3 and r=4 as per user request
            if r_val in [3, 4]:
                plot_df = plot_df[plot_df["s"] <= 15]

            s_vals = sorted(plot_df["s"].unique())
            if not s_vals: continue

            print(f"\n[Plotting Time] Dataset: {ds}, r: {r_val}")
            print(plot_df[plot_df['source'] == 'CBS'][['s', 'source', 'total_sec']].sort_values(['s']).to_string(index=False))

            fig, ax = plt.subplots(figsize=(2, 1.5))
            series_time, first_timeout_idx = {}, {}
            pivot = plot_df.pivot_table(index='s', columns='source', values='total_sec').reindex(s_vals)

            # Identify all 's' values where at least one algorithm succeeded
            # We will use this to mark missing data as failure if others succeeded
            successful_s_set = set()
            for s in s_vals:
                if not pivot.loc[s].isna().all():
                    successful_s_set.add(s)

            for src in _order_series_keys(plot_df["source"].unique()):
                # Even if src is not in pivot columns (completely missing), we might need to mark it as failed
                # But usually it will be in columns if it appeared at least once in the filtered df.
                # If it's not in pivot, it means it has no data for this (ds, r).
                # However, the loop iterates over plot_df["source"].unique(), so it must be in pivot.
                if src not in pivot.columns: continue
                
                ys, first_nan = [], True
                for j, s in enumerate(s_vals):
                    t = pivot.at[s, src]
                    
                    # Check if this algorithm succeeded
                    is_success = pd.notna(t) and t > 0
                    
                    if is_success:
                        ys.append(t)
                    else:
                        # It failed (NaN or 0) OR it was missing but others succeeded
                        # If t is NaN, it's already missing/failed.
                        # We also check if this 's' was solvable by ANY algorithm.
                        # If yes, and this algo failed/missing, we mark it as a timeout point.
                        if first_nan and (s in successful_s_set):
                            ys.append(TIME_YMAX_SEC)
                            first_timeout_idx[src] = j
                            first_nan = False
                        else:
                            ys.append(np.nan)
                series_time[src] = ys

            for lab, ys in series_time.items():
                ax.plot(np.arange(len(s_vals)), ys, label=lab, color=COLOUR_MAP[lab], linestyle=LINESTYLE_MAP[lab],
                        linewidth=2, marker=MARKER_MAP[lab], markersize=MARKER_SIZE_PT,
                        markerfacecolor='none' if HOLLOW_MARKERS else 'auto', markeredgecolor=COLOUR_MAP[lab],
                        markeredgewidth=MARKER_EDGE_WIDTH)

            for lab, j in first_timeout_idx.items():
                ax.scatter(j, TIME_YMAX_SEC * 0.995, s=TIMEOUT_MARKER_SIZE, marker='x', color='red',
                           linewidths=max(1.2, MARKER_EDGE_WIDTH), zorder=7, clip_on=False)

            tick_positions, tick_labels = _get_smart_ticks(s_vals)
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, rotation=60, ha='right')
            ax.set_xlim(-0.5, len(s_vals) - 0.5)
            ax.set_yscale('log')
            ax.set_ylim(bottom=EPS_T, top=TIME_YMAX_SEC)
            ax.set_xlabel("s")
            ax.set_ylabel("Time (s)")
            _style_axes(ax)

            _save_standalone_legend(ax, out_dir, f"{ds}_legend_time", ncol=len(series_time))
            fig.savefig(os.path.join(out_dir, f"{ds}_r{r_val}_time_by_s.eps"), dpi=300, bbox_inches='tight')
            plt.close(fig)


def _plot_memory_compare_by_r(df: pd.DataFrame, out_dir: str):
    """For each dataset and r, draw a line chart for memory vs s across algorithms."""
    if df.empty: return

    # Filter to keep only the requested sources for memory plots
    # Added CBS3 to the list
    target_sources = ['CBS-noHi', 'ARB-noHi', 'Nuclear-NO', 'CBS3']
    
    for ds in sorted(df["dataset_name"].unique()):
        for r_val in sorted(df[df["dataset_name"] == ds]["r"].unique()):
            plot_df = df[(df["dataset_name"] == ds) & (df["r"] == r_val)]

            # Filter for r=3 and r=4 as per user request
            if r_val in [3, 4]:
                plot_df = plot_df[plot_df["s"] <= 15]
            # Filter for r=1 and r=2 as per user request
            elif r_val in [1, 2]:
                plot_df = plot_df[plot_df["s"] <= 30]

            # Apply source filtering
            plot_df = plot_df[plot_df['source'].isin(target_sources)]

            # --- Outlier Filtering ---
            # Remove suspicious low memory values for soc-pokec-relationships
            if ds == 'soc-pokec-relationships':
                plot_df = plot_df[~((plot_df['max_rss_mb'] < 1000))]

            s_vals = sorted(plot_df["s"].unique())
            if not s_vals: continue

            print(f"\n[Plotting Memory] Dataset: {ds}, r: {r_val}")
            # Only print if CBS is present (which it won't be after filtering, so this line is effectively disabled for now)
            # print(plot_df[plot_df['source'] == 'CBS'][['s', 'source', 'max_rss_mb']].sort_values(['s']).to_string(index=False))

            fig, ax = plt.subplots(figsize=(2, 1.5))
            series_mem = {}
            pivot = plot_df.pivot_table(index='s', columns='source', values='max_rss_mb').reindex(s_vals)

            for src in _order_series_keys(plot_df["source"].unique()):
                if src not in pivot.columns: continue
                ys = pivot[src].tolist()
                series_mem[src] = ys

            for lab, ys in series_mem.items():
                # Filter out series that are all NaN
                if all(pd.isna(y) for y in ys): continue

                ax.plot(np.arange(len(s_vals)), ys, label=lab, color=COLOUR_MAP.get(lab, 'black'),
                        linestyle=LINESTYLE_MAP.get(lab, 'solid'), linewidth=2,
                        marker=MARKER_MAP.get(lab, 'o'), markersize=MARKER_SIZE_PT,
                        markerfacecolor='none' if HOLLOW_MARKERS else 'auto',
                        markeredgecolor=COLOUR_MAP.get(lab, 'black'),
                        markeredgewidth=MARKER_EDGE_WIDTH)

            tick_positions, tick_labels = _get_smart_ticks(s_vals)
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, rotation=60, ha='right')
            ax.set_xlim(-0.5, len(s_vals) - 0.5)
            
            # Use symlog for a more "aggressive" log strategy that handles wide ranges better
            ax.set_yscale('symlog', linthresh=100)
            ax.set_ylim(bottom=1) # Set minimum Y value to 1 as requested
            ax.grid(True, which="both", ls="--", alpha=0.3) # Add grid to help read log scale
            
            ax.set_xlabel("s")
            ax.set_ylabel("Memory (MB)")
            _style_axes(ax)

            _save_standalone_legend(ax, out_dir, f"{ds}_legend_mem", ncol=len(series_mem))
            fig.savefig(os.path.join(out_dir, f"{ds}_r{r_val}_memory_by_s.eps"), dpi=300, bbox_inches='tight')
            plt.close(fig)


# --- Data Loading and Processing ---

def _regenerate_csv_from_raw_logs():
    """
    Processes raw experimental logs and saves them as standardized CSV files.
    This is a slow operation and should only be run when raw data changes.
    """
    print("Regenerating data from raw sources (this may be slow)...")

    # As per user request, previously commented-out data sources are now processed.
    # These were likely commented out for performance, not because they are unused.

    # Source for data1.csv
    raw_cbs_main = make_plots.myData(os.path.join(DATA_DIR, "experimentdataNew"))
    _prepare_cbs_data(raw_cbs_main).to_csv(os.path.join(DATA_DIR, "data1.csv"), index=False)

    # Source for data2.csv
    raw_arb = make_plots_bazel.arcData()
    _prepare_arb_data(raw_arb).to_csv(os.path.join(DATA_DIR, "data2.csv"), index=False)
    raw_arb.to_csv(os.path.join(DATA_DIR, "data2_raw.csv"), index=False)

    # Source for data3.csv (Nuclear) - loading mechanism is not clear in this script.
    # The original commented line was:
    # raw_nuclear = arcData(os.path.join(DATA_DIR, "nucleus_experiment.log"))
    # This is likely incorrect as arcData is for ARB. Assuming data3.csv is handled manually.

    # Source for data4.csv (CBS-noHi)
    # UPDATED: Now includes expDataNew.txt as per user request
    raw_cbs_hi_source = make_plots.myData(os.path.join(DATA_DIR, "experimentdataHi"))
    raw_cbs_from_new_txt = make_plots.myData(os.path.join(DATA_DIR, "expDataNew.txt"))
    raw_cbs_hi_combined = pd.concat([raw_cbs_hi_source, raw_cbs_from_new_txt], ignore_index=True)
    
    df_cbs_nohi = _prepare_cbs_data(raw_cbs_hi_combined)
    df_cbs_nohi['source'] = 'CBS-noHi'  # Override source for this dataset
    df_cbs_nohi.to_csv(os.path.join(DATA_DIR, "data4.csv"), index=False)

    # Source for dataAlls.csv
    # UPDATED: Removed expDataNew.txt from here
    raw_cbs_alls = make_plots.myData(os.path.join(DATA_DIR, "experimentdataAlls"))
    _prepare_cbs_data(raw_cbs_alls).to_csv(os.path.join(DATA_DIR, "dataAlls.csv"), index=False)

    # Source for CBS3 (dataCBS3.csv)
    # UPDATED: Added new source for CBS3
    raw_cbs3 = make_plots.myData(os.path.join(DATA_DIR, "expDataone.out"))
    df_cbs3 = _prepare_cbs_data(raw_cbs3)
    df_cbs3['source'] = 'CBS3'
    df_cbs3.to_csv(os.path.join(DATA_DIR, "dataCBS3.csv"), index=False)

    # Other data sources
    raw_scal_new = make_plots.myData(os.path.join(DATA_DIR, "experimentdataScalNew1"))
    _prepare_cbs_data(raw_scal_new).to_csv(os.path.join(DATA_DIR, "dataScalNew.csv"), index=False)

    raw_arb_nohi = make_plots_bazel.arcData(base_root_dir=os.path.join(DATA_DIR, "data_bazel"))
    df_arb_nohi = _prepare_arb_data(raw_arb_nohi)
    df_arb_nohi['source'] = "ARB-noHi"
    df_arb_nohi.to_csv(os.path.join(DATA_DIR, "dataARBnoHi.csv"), index=False)

    print("Finished regenerating and saving data to CSV.")


def _load_processed_data() -> pd.DataFrame:
    """Loads all pre-processed CSV files into a single combined DataFrame."""
    print("Loading pre-processed data from CSV files...")

    def read_csv_safe(filename):
        path = os.path.join(DATA_DIR, filename)
        return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()

    df_cbs = read_csv_safe("data1.csv")
    df_arb = read_csv_safe("data2.csv")
    df_nuclear = read_csv_safe("data3.csv")
    df_cbs_nohi = read_csv_safe("data4.csv")
    df_cbs_alls = read_csv_safe("dataAlls.csv")
    df_arb_nohi = read_csv_safe("dataARBnoHi.csv")
    df_cbs3 = read_csv_safe("dataCBS3.csv") # Load CBS3 data

    df_cbs_alls = df_cbs_alls[~df_cbs_alls["dataset_name"].str.contains(r"\.p", na=False)]
    df_cbs = pd.concat([df_cbs, df_cbs_alls], ignore_index=True)

    # Note: Source assignment logic has been moved to _prepare_cbs_data and
    # _regenerate_csv_from_raw_logs for better separation of concerns.

    combined = pd.concat([df_cbs, df_arb, df_nuclear, df_cbs_nohi, df_arb_nohi, df_cbs3], ignore_index=True)
    for col in ['r', 's']: combined[col] = pd.to_numeric(combined[col], errors='coerce').astype('Int64')
    return _dedup_min_time_per_source(combined)


def _apply_general_filters(df: pd.DataFrame) -> pd.DataFrame:
    """Applies general filters to the combined dataframe before plotting."""
    if df.empty: return df

    s_bounds = {'com-dblp': 114, 'com-youtube': 17, 'web-Google': 44, 'web-Stanford': 61,
                'soc-pokec-relationships': 29}
    for ds, s_max in s_bounds.items():
        df = df[~((df['dataset_name'] == ds) & (df['s'] > s_max))]

    df = df[(df['total_sec'] <= TIME_YMAX_SEC) & (df['r'] != df['s'])]
    return df


def main():
    """Main script for generating comparison plots."""
    if REGENERATE_FROM_RAW:
        _regenerate_csv_from_raw_logs()

    combined_df = _load_processed_data()
    combined_df = _apply_general_filters(combined_df)

    # --- Analysis ---
    _analyze_cbs_vs_nohi_speedup(combined_df.copy())

    # --- Plotting ---
    _plot_cbs_vs_nohi_long_comparison(combined_df.copy(), COMPARE_OUT)

    # --- Experimental Filters (preserved from original script) ---
    # This section contains various filters that were commented out in the original script.
    # They are likely used for experimenting with different views of the data.
    # exit(0)

    # Example: Keep only 'noHi' or 'NO' variants and rename them
    # combined_df = combined_df[combined_df['source'].isin(['CBS-noHi', 'ARB-noHi', 'Nuclear-NO'])]
    # combined_df['source'] = combined_df['source'].map({'CBS-noHi': 'CBS', 'ARB-noHi': 'ARB', 'Nuclear-NO': 'Nuclear'})

    # Example: Filter by r or s values
    # combined_df = combined_df[(combined_df['r'] <= 2)]
    # combined_df = combined_df[(combined_df['r'] > 2) & (combined_df['s'] <= 15)]

    # --- Final Plots and Summary ---
    print("\nAlgorithms present in the final filtered dataset:", sorted(combined_df['source'].unique()))
    _plot_compare_by_r(combined_df, COMPARE_OUT)
    _plot_memory_compare_by_r(combined_df, COMPARE_OUT)
    # _print_speedups_by_dataset(combined, DATA_DIR)
    print("\nDone.")


if __name__ == "__main__":
    os.makedirs(COMPARE_OUT, exist_ok=True)
    main()
