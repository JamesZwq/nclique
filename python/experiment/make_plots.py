#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pathlib
from typing import List, Dict, Any
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

MAX_HEAP_PLOTS_PER_R = 30  # 每个 (dataset, r_small) 只取 ~10 个 S 采样来画 heap 图

RAW_PATH = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdataAlls"
OUT_DIR  = pathlib.Path("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/image/")

# -------------------------
# Helpers
# -------------------------

def ds_name(p: str):
    if not isinstance(p, str) or not p:
        return p
    base = os.path.basename(p.strip())
    # Repeatedly strip known extensions (handles things like ".mtx.gz")
    known_exts = {".edges", ".adj", ".txt", ".mtx", ".csv", ".tsv", ".bin", ".gz"}
    while True:
        root, ext = os.path.splitext(base)
        if ext and ext.lower() in known_exts and root:
            base = root
            continue
        break
    return base

def is_run_start(line: str) -> bool:
    # e.g. "================== RUN =================="
    return bool(re.match(r"^=+\s*RUN\s*=+\s*$", line.strip(), re.I))

def is_run_end(line: str) -> bool:
    # e.g. "================ END OF RUN =================="
    return bool(re.match(r"^=+\s*END OF RUN\s*=+\s*$", line.strip(), re.I))

def split_runs(text: str) -> List[str]:
    """
    Robustly split the log into run blocks using a regex pattern that captures
    content between a RUN marker and the next END OF RUN marker.
    """
    # This pattern finds all non-overlapping occurrences of blocks.
    # - It starts matching at a RUN line.
    # - (.*?) non-greedily captures all characters (including newlines due to re.DOTALL).
    # - It stops at the first END OF RUN line it finds.
    pattern = re.compile(
        r"^=+\s*RUN\s*=+\s*\n(.*?)\n^=+\s*END OF RUN\s*=+",
        re.DOTALL | re.MULTILINE | re.IGNORECASE
    )
    return pattern.findall(text)

# Patterns for fields & metrics (tolerant)
FIELD_PAT = {
    "timestamp": re.compile(r"^timestamp\s*:\s*(.+)$", re.M | re.I),
    "dataset":   re.compile(r"^dataset\s*:\s*(.+)$", re.M | re.I),
    "r":         re.compile(r"^s\s*:\s*([0-9]+)\s*$", re.M | re.I),
    "s":         re.compile(r"^r\s*:\s*([0-9]+)\s*$", re.M | re.I),
    "cmd":       re.compile(r"^cmd\s*:\s*(.+)$", re.M | re.I),
}

METRIC_PAT = {
    # Tree Build took: 127918 ms
    "tree_build_ms":   re.compile(r"Tree\s*Build\s*took:\s*([\d.]+)\s*ms", re.I),
    # clique Index build took: 4263.32 ms (sometimes "Clique index")
    "clique_index_ms": re.compile(r"clique\s*index\s*build\s*took:\s*([\d.]+)\s*ms", re.I),
    # Accept "time: 58055481 ms" or "NucleusCoreDecomposition took: 5.80907e+07 ms"
    "prog_time_ms":    re.compile(r"(?:\btime\s*:\s*|NucleusCoreDecomposition\s+took\s*:\s*)([\d.e+]+)\s*ms\b", re.I),
    # nun Leaf: 907668 / num Leaf / nu Leaf
    "num_leaves":      re.compile(r"\b(?:nu[mn]?|num)\s*Leaf\s*:\s*([0-9]+)", re.I),
    # Maximum resident set size (kbytes): 3895968
    "max_rss_kb":      re.compile(r"Maximum\s+resident\s+set\s+size\s*\(kbytes\)\s*:\s*([0-9]+)", re.I),
    # Exit status: 0
    "exit_status":     re.compile(r"Exit\s+status\s*:\s*([0-9]+)", re.I),
}

# Lines like: "minCore: 0.00, heap size: 5402"
HEAP_LINE_PAT = re.compile(r"minCore\s*:\s*([0-9.]+)\s*,\s*heap\s*size\s*:\s*([0-9]+)", re.I)

def parse_block(block: str) -> Dict[str, Any]:
    row: Dict[str, Any] = {k: None for k in ["dataset","s","r","cmd"]}
    for k, pat in FIELD_PAT.items():
        m = pat.search(block)
        if m:
            val = m.group(1).strip()
            if k in ("s","r"):
                try:
                    row[k] = int(val)
                except Exception:
                    row[k] = None
            else:
                row[k] = val

    for k, pat in METRIC_PAT.items():
        m = pat.search(block)
        if not m:
            row[k] = None
            continue
        val = m.group(1)
        if k in ("num_leaves","max_rss_kb","exit_status"):
            try:
                row[k] = int(val)
            except Exception:
                row[k] = None
        else:
            try:
                row[k] = float(val)
            except Exception:
                row[k] = None

    if re.search(r"Command\s+terminated\s+by\s+signal", block, flags=re.I):
        row["exit_status"] = -1
        # print(f"[WARN] Command terminated by signal; setting exit_status=-1 for dataset: {row.get('dataset')} cmd: {row.get('cmd')}")

    # Parse heap-series lines (if any)
    heap_series = []
    for mc_str, hs_str in HEAP_LINE_PAT.findall(block):
        try:
            heap_series.append((float(mc_str), int(hs_str)))
        except Exception:
            continue
    row["heap_series"] = heap_series if heap_series else None

    return row

def myData(path: str) -> pd.DataFrame:
    # Read the entire file once
    text = pathlib.Path(path).read_text(encoding="utf-8", errors="replace")
    
    # Get all potential blocks that match the RUN...END pattern
    potential_blocks = split_runs(text)
    
    # Filter for compliant blocks based on the new strict criteria
    compliant_blocks = []
    for block in potential_blocks:
        begin_marker = '=========================begin========================='
        time_marker = '---------------- /bin/time -v ----------------'
        
        # Check for exactly one occurrence of each marker
        if block.count(begin_marker) != 1 or block.count(time_marker) != 1:
            continue
            
        # Check for the correct order
        try:
            begin_pos = block.index(begin_marker)
            time_pos = block.index(time_marker, begin_pos) # Search for time_marker after begin_marker
            compliant_blocks.append(block)
        except ValueError:
            # This block doesn't have the markers in the correct order, so we skip it.
            continue
            
    # Parse only the compliant blocks
    rows = [parse_block(b) for b in compliant_blocks]

    # Create DataFrame and ensure we only have successful runs
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df[df['exit_status'] == 0]

    # Basic cleaning
    if "dataset" in df.columns:
        df["dataset_name"] = df["dataset"].map(ds_name)
    else:
        df["dataset_name"] = None

    # normalize numeric columns
    for c in ("s","r","tree_build_ms","clique_index_ms","prog_time_ms"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    return df

def dedupe_best(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # Only keep successful runs if exit_status exists
    if "exit_status" in df.columns:
        ok = df["exit_status"].isna() | (df["exit_status"] == 0)
        df = df[ok]

    # Keep rows that have at least two of the three timing metrics present.
    # (Drop rows with 2 or more NULLs among tree_build_ms / clique_index_ms / prog_time_ms)
    required = ["tree_build_ms","clique_index_ms","prog_time_ms"]
    have = [c for c in required if c in df.columns]
    if have:
        nonnull_cnt = df[have].notna().sum(axis=1)
        df = df[nonnull_cnt >= 2].copy()
        # For plotting/stacking, treat the at-most-one missing metric as 0
        for c in have:
            df[c] = df[c].fillna(0.0)

    # Deduplicate by (dataset, small-parameter, big-parameter): prefer most-complete row
    def score(row):
        return int(row["tree_build_ms"] > 0) + int(row["clique_index_ms"] > 0) + int(row["prog_time_ms"] > 0)

    # Identify the parameter columns after normalization:
    # - small parameter in {1,2,3}: "r_small" (fallback to "r")
    # - large parameter in [3,30]: "S" (fallback to "s")
    key_small = "r_small" if "r_small" in df.columns else ("r" if "r" in df.columns else None)
    key_big   = "S"       if "S"       in df.columns else ("s" if "s" in df.columns else None)

    # Drop rows missing key identifiers
    if key_small is None or key_big is None:
        # If parameters cannot be identified, at least keep dataset name
        df = df.dropna(subset=["dataset_name"])
        return df
    df = df.dropna(subset=["dataset_name", key_small, key_big])

    if len(df):
        df["_score"] = df.apply(score, axis=1)
        # Keep the row with highest score per (dataset, small, big)
        df = (
            df.sort_values(["dataset_name", key_small, key_big, "_score"],
                           ascending=[True, True, True, False])
              .drop_duplicates(subset=["dataset_name", key_small, key_big], keep="first")
              .drop(columns=["_score"])
        )
    return df

def normalize_params(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of df with two canonical columns:
    - 'r_small' for the small clique parameter in {1,2,3}
    - 'S' for the large clique size in [3,30]
    The raw logs sometimes swap meanings of the printed fields 's' and 'r';
    this normalizer auto-detects and renames accordingly.
    """
    df = df.copy()
    # Coerce numeric
    for c in ("s", "r"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Heuristic: whichever column has max <= 3 is the small parameter
    s_col, S_col = "s", "r"
    s_max = pd.to_numeric(df.get("s"), errors="coerce").max()
    r_max = pd.to_numeric(df.get("r"), errors="coerce").max()
    if pd.notna(s_max) and pd.notna(r_max):
        if s_max <= 4 and r_max >= 4:
            s_col, S_col = "s", "r"
        elif r_max <= 4 and s_max >= 4:
            s_col, S_col = "r", "s"
        else:
            # Fallback: keep default, but the downstream filter will likely drop nonconforming rows
            s_col, S_col = "s", "r"
    df = df.rename(columns={s_col: "r_small", S_col: "S"})
    return df

def plot_one(ds_name_str, r_val, g, out_dir: pathlib.Path):
    """Plot stacked bars over S (large clique size) for a fixed small r=r_val."""
    # Keep only S in [3, 30]
    if "S" not in g.columns:
        return None
    g = g[(g["S"] >= 3) & (g["S"] <= 30)].copy()
    if g.empty:
        return None
    g = g.sort_values("S")
    s_vals = g["S"].astype(int).tolist()

    tree = g.get("tree_build_ms", pd.Series([0.0] * len(g))).values
    cli  = g.get("clique_index_ms", pd.Series([0.0] * len(g))).values
    run  = g.get("prog_time_ms", pd.Series([0.0] * len(g))).values

    width = max(6.0, 0.45 * len(s_vals) + 2.0)
    plt.figure(figsize=(width, 4.6))
    x = np.arange(len(s_vals))

    # Stacked bars: Tree (bottom) + CliqueIndex (middle, if present) + Running (top)
    plt.bar(x, tree, label="Tree build (ms)")
    bottom2 = tree
    if np.any(np.nan_to_num(cli) > 0):
        plt.bar(x, cli, bottom=tree, label="Clique index (ms)")
        bottom2 = tree + cli
    plt.bar(x, run, bottom=bottom2, label="Running time (ms)")

    plt.xticks(x, s_vals, rotation=0)
    plt.xlabel("S")
    plt.ylabel("Time (ms)")
    plt.title(f"{ds_name_str} — r={int(r_val)}")
    plt.legend()
    plt.tight_layout()

    out_path = out_dir / f"{ds_name_str}_r{int(r_val)}.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    return out_path
def plot_mem_one(ds_name_str, r_val, g, out_dir: pathlib.Path):
    """Plot bars of peak memory (Max RSS) over S for a fixed r=r_val."""
    if "S" not in g.columns or "max_rss_kb" not in g.columns:
        return None
    g = g[(g["S"] >= 3) & (g["S"] <= 30)].copy()
    g = g.dropna(subset=["max_rss_kb"])
    if g.empty:
        return None

    g = g.sort_values("S")
    s_vals = g["S"].astype(int).tolist()

    # Convert kbytes -> GiB
    mem_gib = (g["max_rss_kb"].astype(float) / (1024.0 * 1024.0)).values

    width = max(6.0, 0.45 * len(s_vals) + 2.0)
    plt.figure(figsize=(width, 4.6))
    x = np.arange(len(s_vals))

    # Single-segment bars: peak process memory
    plt.bar(x, mem_gib, label="Max RSS (GiB)")

    plt.xticks(x, s_vals, rotation=0)
    plt.xlabel("S")
    plt.ylabel("Memory (GiB)")
    plt.title(f"{ds_name_str} — r={int(r_val)} (memory)")
    plt.legend()
    plt.tight_layout()

    out_path = out_dir / f"{ds_name_str}_r{int(r_val)}_mem.png"
    plt.savefig(out_path, dpi=150)
    plt.close()
    return out_path

def plot_heap_series(ds_name_str, r_val, S_val, series, out_dir: pathlib.Path):
    """Plot per-step heap sizes as a bar chart for one (dataset, r, S). Also writes a CSV."""
    if not series:
        return None
    # 忽视 heap_size == 0 的点（logscale 下为无效），只保留 >0
    filtered = [(mc, h) for mc, h in series if (h is not None and h > 0)]
    if not filtered:
        return None

    steps = np.arange(len(filtered))
    mincores = [mc for mc, _ in filtered]
    heaps = [h for _, h in filtered]

    width = max(6.0, 0.35 * len(filtered) + 2.0)
    plt.figure(figsize=(width, 4.6))

    plt.bar(steps, heaps, label="Heap size")
    plt.yscale("log")

    # Keep x ticks readable: ~<=20 ticks
    try:
        import math
        step = max(1, int(math.ceil(len(filtered) / 20.0)))
    except Exception:
        step = max(1, int((len(filtered) + 19) // 20))
    xticks = steps[::step]
    plt.xticks(xticks, [str(int(i)) for i in xticks], rotation=0)

    plt.xlabel("Step")
    plt.ylabel("Heap size (log)")
    plt.title(f"{ds_name_str} — r={int(r_val)}, S={int(S_val)} (heap)")
    plt.legend()
    plt.tight_layout()

    out_path = out_dir / f"{ds_name_str}_r{int(r_val)}_S{int(S_val)}_heap.png"
    plt.savefig(out_path, dpi=150)
    plt.close()

    # Also write a CSV alongside the plot for downstream analysis
    csv_path = out_dir / f"{ds_name_str}_r{int(r_val)}_S{int(S_val)}_heap.csv"
    pd.DataFrame({"step": steps, "minCore": mincores, "heap_size": heaps}).to_csv(csv_path, index=False)
    print(f"Wrote heap series CSV: {csv_path}")
    return out_path

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Parse and clean
    df = myData(RAW_PATH)
    df1 = myData("/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdata")
    if df.empty:
        print(f"[WARN] No runs parsed from: {RAW_PATH}")
        return

    df = pd.concat([df, df1], ignore_index=True)

    # Canonicalize parameter columns: 'r_small' in {1,2,3}, 'S' in [3,30]
    df = normalize_params(df)
    print(df.columns.tolist())
    # Keep only the three small-clique settings we care about
    if "r_small" not in df.columns:
        print("[WARN] Could not infer small-clique parameter from logs; no plots generated.")
        return
    df = df[df["r_small"].isin([1, 2, 3, 4])]
    if df.empty:
        print("[WARN] No entries with r_small in {1,2,3} in the experiment data.")
        return

    df = dedupe_best(df)

    generated = []
    # Group by dataset then r_small
    for ds, gds in df.groupby("dataset_name", dropna=True):
        for r_val in (1, 2, 3, 4):
            gr = gds[gds["r_small"] == r_val]
            if gr.empty:
                continue
            p = plot_one(ds, r_val, gr, OUT_DIR)
            if p:
                generated.append(str(p))

    # Write an index file for convenience
    idx = OUT_DIR / "_plots_index.csv"
    pd.Series(generated, name="plot_path").to_csv(idx, index=False)
    print(f"Saved {len(generated)} plots under: {OUT_DIR}")
    print(f"Index: {idx}")
    # ----- Memory plots (Max RSS) -----
    generated_mem = []
    for ds, gds in df.groupby("dataset_name", dropna=True):
        for r_val in (1, 2, 3, 4):
            gr = gds[gds["r_small"] == r_val]
            if gr.empty:
                continue
            p = plot_mem_one(ds, r_val, gr, OUT_DIR)
            if p:
                generated_mem.append(str(p))

    idx_mem = OUT_DIR / "_plots_mem_index.csv"
    pd.Series(generated_mem, name="plot_path").to_csv(idx_mem, index=False)
    print(f"Saved {len(generated_mem)} memory plots under: {OUT_DIR}")
    print(f"Memory index: {idx_mem}")

    # ----- Heap-series bar plots (sample ~10 S per (dataset, r_small)) -----
    generated_heap = []
    # 只保留含 heap_series 与 S 的行
    df_heap = df.copy()
    if "heap_series" not in df_heap.columns:
        df_heap["heap_series"] = None
    df_heap = df_heap[~df_heap["heap_series"].isna()]
    if "S" in df_heap.columns:
        df_heap = df_heap[~df_heap["S"].isna()]
    else:
        df_heap = pd.DataFrame(columns=df_heap.columns)  # 无 S 列则直接空
    if not df_heap.empty:
        # 分组：每个 (dataset, r_small) 里按 S 排序，然后等距采样最多 MAX_HEAP_PLOTS_PER_R 个
        for (ds, r_val), ggrp in df_heap.groupby(["dataset_name", "r_small"], dropna=True):
            ggrp = ggrp.sort_values("S")
            if ggrp.empty:
                continue
            n = int(min(MAX_HEAP_PLOTS_PER_R, len(ggrp)))
            if n <= 0:
                continue
            # 等距索引（去重并保持有序）
            idxs = np.linspace(0, len(ggrp) - 1, num=n)
            idxs = np.round(idxs).astype(int)
            idxs = sorted(set(int(i) for i in idxs))
            gsel = ggrp.iloc[idxs]
            for _, row in gsel.iterrows():
                series = row.get("heap_series", None)
                S_val = row.get("S", None)
                if not series or pd.isna(S_val):
                    continue
                p = plot_heap_series(str(ds), int(r_val), int(S_val), series, OUT_DIR)
                if p:
                    generated_heap.append(str(p))

    idx_heap = OUT_DIR / "_plots_heap_index.csv"
    pd.Series(generated_heap, name="plot_path").to_csv(idx_heap, index=False)
    print(f"Saved {len(generated_heap)} heap plots under: {OUT_DIR}")
    print(f"Heap index: {idx_heap}")
    # Simple summary to stdout for quick sanity check
    if len(generated) == 0:
        print("[HINT] If you expected plots, check that your experimentdata contains r ∈ {1,2,3}.\n"
              "       Try grepping for lines like 'r         : 1' etc.")

if __name__ == "__main__":
    main()
