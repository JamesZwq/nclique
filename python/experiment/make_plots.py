#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import pathlib
from typing import List, Dict, Any
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

RAW_PATH = "/Users/zhangwenqian/UNSW/pivoter/python/experiment/data/experimentdata"
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
    """Robustly split the log into run blocks between RUN/END markers.
    If END is missing for the last block, take till EOF."""
    lines = text.splitlines(True)  # keep line breaks
    blocks = []
    cur = []
    inside = False
    for ln in lines:
        if is_run_start(ln):
            # flush previous if somehow unclosed
            if inside and cur:
                blocks.append("".join(cur))
                cur = []
            inside = True
            cur = []  # start fresh (exclude the header line itself)
            continue
        if is_run_end(ln):
            if inside:
                blocks.append("".join(cur))
            inside = False
            cur = []
            continue
        if inside:
            cur.append(ln)
    if inside and cur:
        # last block without explicit END
        blocks.append("".join(cur))
    return blocks

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
    # Accept either "time: 1912 ms" or "NucleusCoreDecomposition took: 1912 ms"
    "prog_time_ms":    re.compile(r"(?:\btime\s*:\s*|NucleusCoreDecomposition\s+took\s*:\s*)([\d.]+)\s*ms\b", re.I),
    # nun Leaf: 907668 / num Leaf / nu Leaf
    "num_leaves":      re.compile(r"\b(?:nu[mn]?|num)\s*Leaf\s*:\s*([0-9]+)", re.I),
    # Maximum resident set size (kbytes): 3895968
    "max_rss_kb":      re.compile(r"Maximum\s+resident\s+set\s+size\s*\(kbytes\)\s*:\s*([0-9]+)", re.I),
    # Exit status: 0
    "exit_status":     re.compile(r"Exit\s+status\s*:\s*([0-9]+)", re.I),
}

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
    return row

def myData(path: str) -> pd.DataFrame:
    # Read once
    text = pathlib.Path(path).read_text(encoding="utf-8", errors="replace")
    blocks = split_runs(text)
    rows = [parse_block(b) for b in blocks]
    # remove if exit_status is not 0
    rows = [r for r in rows if r.get("exit_status") == 0]
    df = pd.DataFrame(rows)

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
        if s_max <= 3 and r_max >= 3:
            s_col, S_col = "s", "r"
        elif r_max <= 3 and s_max >= 3:
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
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Parse and clean
    df = myData(RAW_PATH)
    if df.empty:
        print(f"[WARN] No runs parsed from: {RAW_PATH}")
        return

    # Canonicalize parameter columns: 'r_small' in {1,2,3}, 'S' in [3,30]
    df = normalize_params(df)
    print(df.columns.tolist())
    # Keep only the three small-clique settings we care about
    if "r_small" not in df.columns:
        print("[WARN] Could not infer small-clique parameter from logs; no plots generated.")
        return
    df = df[df["r_small"].isin([1, 2, 3])]
    if df.empty:
        print("[WARN] No entries with r_small in {1,2,3} in the experiment data.")
        return

    df = dedupe_best(df)

    generated = []
    # Group by dataset then r_small
    for ds, gds in df.groupby("dataset_name", dropna=True):
        for r_val in (1, 2, 3):
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
        for r_val in (1, 2, 3):
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
    # Simple summary to stdout for quick sanity check
    if len(generated) == 0:
        print("[HINT] If you expected plots, check that your experimentdata contains r ∈ {1,2,3}.\n"
              "       Try grepping for lines like 'r         : 1' etc.")

if __name__ == "__main__":
    main()