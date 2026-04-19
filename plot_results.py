#!/usr/bin/env python3
"""
Fetch results from tods1/tods2, merge, and plot heatmaps.
Usage: python3 plot_results.py
"""

import subprocess, csv, sys, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from collections import defaultdict
from pathlib import Path

SERVERS = ["tods1", "tods2"]
LOCAL_CSV = "bench_v3_all_results.csv"
REMOTE_CSV = "~/nclique/bench_v3_all_results.csv"

ALGOS = ["REF", "ST", "V3Fast", "V3LM", "V3Fast_NP", "V3Fast_NoCPI", "V3H", "V3HC", "V3", "V3_NP"]
ALGO_TITLES = {
    "REF": "REF (baseline)",
    "ST": "ST (CPI peeling)",
    "V3Fast": "V3Fast (optimized + Private Cloud)",
    "V3LM": "V3LM (V3Fast, low-memory engineering)",
    "V3Fast_NP": "V3Fast (no Private Cloud) — PC ablation",
    "V3Fast_NoCPI": "V3Fast (no CPI, enumerate) — CPI ablation",
    "V3H": "V3H (V3Fast + tuple-based hierarchy)",
    "V3HC": "V3HC (V3Fast + class-based hierarchy)",
    "V3": "V3 (legacy, pre-optimization)",
    "V3_NP": "V3_NP (legacy, pre-optimization)",
}
TIMEOUT_VAL = -0.5
OOM_VAL = -1.5

# ============ Fetch & Merge ============
def fetch():
    rows = []
    header = None
    for server in SERVERS:
        print(f"  Fetching from {server}...", end=" ", flush=True)
        try:
            r = subprocess.run(
                ["ssh", server, f"cat {REMOTE_CSV}"],
                capture_output=True, text=True, timeout=15
            )
            if r.returncode != 0:
                print("FAIL")
                continue
            lines = r.stdout.strip().split("\n")
            reader = csv.DictReader(lines)
            if header is None:
                header = reader.fieldnames
            server_rows = list(reader)
            print(f"{len(server_rows)} rows")
            rows.extend(server_rows)
        except Exception as e:
            print(f"ERROR: {e}")

    # Also load local CSV if exists
    if os.path.exists(LOCAL_CSV):
        local_rows = list(csv.DictReader(open(LOCAL_CSV)))
        print(f"  Local: {len(local_rows)} rows")
        rows.extend(local_rows)

    # Dedup by (graph, r, s, algo) — keep last
    seen = {}
    for row in rows:
        key = (row["graph"], row["r"], row["s"], row["algo"])
        seen[key] = row
    rows = list(seen.values())
    print(f"  Total: {len(rows)} unique results")
    return rows, header

# ============ Parse ============
def parse(rows):
    data = {}
    for row in rows:
        try:
            g, r, s, algo = row["graph"], int(row["r"]), int(row["s"]), row["algo"]
            st = row["status"]
            wall = float(row.get("wall_ms") or -1)
            mem = float(row.get("mem_kB") or -1)
            if st == "OK" and wall >= 0:
                data[(g, algo, r, s)] = (wall, mem, "OK")
            elif st in ("OOM",):
                data[(g, algo, r, s)] = (-1, -1, "OOM")
            elif st in ("TIMEOUT", "SKIP_TIMEOUT") or st.startswith("ERROR"):
                data[(g, algo, r, s)] = (-1, -1, "TIMEOUT")
        except:
            pass
    return data

# ============ Summary ============
def print_summary(rows):
    by_ga = defaultdict(lambda: defaultdict(int))
    for row in rows:
        by_ga[(row["graph"], row["algo"])][row["status"]] += 1

    # Load max clique sizes from server cache files
    max_clique = defaultdict(int)
    for server in SERVERS:
        try:
            r = subprocess.run(
                ["ssh", server, "cat ~/nclique/bench_v3_max_cliques.json"],
                capture_output=True, text=True, timeout=10
            )
            if r.returncode == 0:
                import json
                for g, mc in json.loads(r.stdout).items():
                    max_clique[g] = max(max_clique[g], mc)
        except: pass
    # Fallback: use max s from data if cache unavailable
    if not max_clique:
        for row in rows:
            try:
                g, s = row["graph"], int(row["s"])
                max_clique[g] = max(max_clique[g], s)
            except: pass

    # Progress per graph
    print(f"\n{'Graph':>15} {'MaxClq':>6}  {'Total':>7}  Progress")
    print("-" * 55)
    for g in sorted(set(row["graph"] for row in rows)):
        mc = max_clique[g]
        total_jobs = sum(s - 3 for s in range(4, mc + 1)) * len(set(row["algo"] for row in rows if row["graph"] == g))
        done_jobs = sum(1 for row in rows if row["graph"] == g)
        pct = 100.0 * done_jobs / total_jobs if total_jobs > 0 else 0
        bar_len = 20
        filled = int(bar_len * pct / 100)
        bar = "█" * filled + "░" * (bar_len - filled)
        print(f"{g:>15} {mc:>6}  {done_jobs:>5}/{total_jobs:<5}  {bar} {pct:.1f}%")

    # Per algo breakdown
    print(f"\n{'Graph':>15} {'Algo':>6} {'OK':>6} {'TO':>5} {'OOM':>4} {'ERR':>4} {'SKIP':>5}")
    print("-" * 52)
    for (g, a) in sorted(by_ga):
        s = by_ga[(g, a)]
        ok = s.get("OK", 0)
        to = s.get("TIMEOUT", 0)
        oom = s.get("OOM", 0)
        sk = s.get("SKIP_TIMEOUT", 0)
        err = sum(v for k, v in s.items() if k.startswith("ERROR"))
        if ok + to + oom + err + sk > 0:
            print(f"{g:>15} {a:>6} {ok:>6} {to:>5} {oom:>4} {err:>4} {sk:>5}")

# ============ Plot ============
def plot_heatmaps(data):
    graphs = sorted(set(g for g, a, r, s in data))

    for graph in graphs:
        algos_with_data = [a for a in ALGOS if any(k[0] == graph and k[1] == a for k in data)]
        if not algos_with_data:
            continue

        n_algos = len(algos_with_data)
        fig, axes = plt.subplots(2, n_algos, figsize=(5.5 * n_algos, 11))
        if n_algos == 1:
            axes = axes.reshape(2, 1)
        fig.suptitle(f"{graph}: Wall Time & Memory\n(black = timeout/error, purple = OOM, white = not run)",
                     fontsize=14, fontweight='bold')

        metrics = [
            (0, 'RdYlGn_r', 2.5, 6.5, [3, 4, 5, 6], ['1s', '10s', '100s', '1000s'], 'Wall Time'),
            (1, 'YlOrRd', 3, 8, [3, 4, 5, 6, 7, 8], ['1MB', '10MB', '100MB', '1GB', '10GB', '100GB'], 'Memory'),
        ]

        for ai, algo in enumerate(algos_with_data):
            r_list = sorted(set(r for g, a, r, s in data if g == graph and a == algo))
            s_list = sorted(set(s for g, a, r, s in data if g == graph and a == algo))

            for ri, (midx, cmap_name, vmin, vmax, ticks, tlabels, suffix) in enumerate(metrics):
                ax = axes[ri, ai]
                if not r_list or not s_list:
                    ax.set_title(f"{ALGO_TITLES.get(algo, algo)}\n{suffix} (no data)")
                    ax.axis('off')
                    continue

                mat = np.full((len(r_list), len(s_list)), np.nan)
                for i, rv in enumerate(r_list):
                    for j, sv in enumerate(s_list):
                        if rv >= sv:
                            continue
                        val = data.get((graph, algo, rv, sv))
                        if val is None:
                            continue
                        if val[2] == "TIMEOUT":
                            mat[i, j] = TIMEOUT_VAL
                        elif val[2] == "OOM":
                            mat[i, j] = OOM_VAL
                        elif midx == 0 and val[0] > 0:
                            mat[i, j] = np.log10(val[0])
                        elif midx == 1 and val[1] > 0:
                            mat[i, j] = np.log10(val[1])

                cmap = plt.colormaps[cmap_name].resampled(256)
                im = ax.imshow(mat, aspect='auto', origin='lower', cmap=cmap,
                               vmin=vmin, vmax=vmax)

                # Black overlay for TIMEOUT
                to_mask = (mat == TIMEOUT_VAL)
                if to_mask.any():
                    overlay = np.ma.masked_where(~to_mask, np.ones_like(mat))
                    ax.imshow(overlay, aspect='auto', origin='lower',
                              cmap=mcolors.ListedColormap(['#1a1a1a']),
                              vmin=0, vmax=1, alpha=0.95)
                # Purple overlay for OOM
                oom_mask = (mat == OOM_VAL)
                if oom_mask.any():
                    overlay = np.ma.masked_where(~oom_mask, np.ones_like(mat))
                    ax.imshow(overlay, aspect='auto', origin='lower',
                              cmap=mcolors.ListedColormap(['#8B008B']),
                              vmin=0, vmax=1, alpha=0.95)

                rs = max(1, len(r_list) // 10)
                ss = max(1, len(s_list) // 10)
                ax.set_yticks(range(0, len(r_list), rs))
                ax.set_yticklabels([r_list[i] for i in range(0, len(r_list), rs)])
                ax.set_xticks(range(0, len(s_list), ss))
                ax.set_xticklabels([s_list[i] for i in range(0, len(s_list), ss)])
                ax.set_xlabel("s")
                ax.set_ylabel("r")
                ax.set_title(f"{ALGO_TITLES.get(algo, algo)}\n{suffix}",
                             fontsize=10, fontweight='bold')

                if ai == n_algos - 1:
                    cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
                    cbar.set_ticks(ticks)
                    cbar.set_ticklabels(tlabels)

        plt.tight_layout(rect=[0, 0, 0.95, 0.93])
        outf = f"bench_v3_{graph}_heatmap.png"
        plt.savefig(outf, dpi=150, bbox_inches='tight')
        print(f"  Saved {outf}")
        plt.close()

# ============ Main ============
def main():
    print("=== Fetching results ===")
    rows, header = fetch()
    if not rows:
        print("No results found!")
        return

    print_summary(rows)

    print("\n=== Plotting ===")
    data = parse(rows)
    plot_heatmaps(data)
    print("\nDone.")

if __name__ == "__main__":
    main()