#!/usr/bin/env python3
"""Plot V3 benchmark results: heatmaps and speedup charts."""

import csv, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

CSV = "/tmp/bench_v3_results.csv"

rows = list(csv.DictReader(open(CSV)))

# Parse wall_ms for OK results
data = {}  # (graph, algo, r, s) -> wall_ms
for row in rows:
    g, r, s, algo, status = row["graph"], int(row["r"]), int(row["s"]), row["algo"], row["status"]
    wall = row.get("wall_ms", "")
    if status == "OK" and wall:
        data[(g, algo, r, s)] = float(wall)
    elif status in ("TIMEOUT", "OOM", "SKIP_TIMEOUT"):
        data[(g, algo, r, s)] = -1  # mark as failed

graphs_with_data = sorted(set(g for g, _, _, _ in data))

for graph in graphs_with_data:
    # Collect all (r,s) for this graph
    all_rs = sorted(set((r, s) for g, a, r, s in data if g == graph))
    if not all_rs:
        continue

    r_vals = sorted(set(r for r, s in all_rs))
    s_vals = sorted(set(s for r, s in all_rs))

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    fig.suptitle(f"{graph} — Wall Clock Time (ms)", fontsize=14, fontweight='bold')

    for ax_idx, algo in enumerate(["ST", "V3", "V3_NP"]):
        ax = axes[ax_idx]
        # Build heatmap matrix
        mat = np.full((len(r_vals), len(s_vals)), np.nan)
        for i, r in enumerate(r_vals):
            for j, s in enumerate(s_vals):
                if r >= s:
                    continue
                val = data.get((graph, algo, r, s))
                if val is not None and val > 0:
                    mat[i, j] = np.log10(val)
                elif val == -1:
                    mat[i, j] = np.nan  # timeout/oom

        im = ax.imshow(mat, aspect='auto', origin='lower', cmap='RdYlGn_r',
                       vmin=0, vmax=6)  # 1ms to 1000s
        ax.set_xlabel("s")
        ax.set_ylabel("r")
        ax.set_title(algo)

        # Set tick labels (show subset to avoid crowding)
        r_step = max(1, len(r_vals) // 10)
        s_step = max(1, len(s_vals) // 10)
        ax.set_yticks(range(0, len(r_vals), r_step))
        ax.set_yticklabels([r_vals[i] for i in range(0, len(r_vals), r_step)])
        ax.set_xticks(range(0, len(s_vals), s_step))
        ax.set_xticklabels([s_vals[i] for i in range(0, len(s_vals), s_step)])

    cbar = fig.colorbar(im, ax=axes, shrink=0.8, label='log10(wall_ms)')
    cbar.set_ticks([0, 1, 2, 3, 4, 5, 6])
    cbar.set_ticklabels(['1ms', '10ms', '100ms', '1s', '10s', '100s', '1000s'])

    plt.tight_layout()
    outf = f"bench_v3_{graph}_heatmap.png"
    plt.savefig(outf, dpi=150, bbox_inches='tight')
    print(f"  Saved {outf}")
    plt.close()

    # --- Speedup plot: V3/ST for each (r,s) ---
    fig2, ax2 = plt.subplots(figsize=(12, 6))
    fig2.suptitle(f"{graph} — V3 Speedup over ST (wall clock)", fontsize=14)

    for r in r_vals:
        ss = []
        speedups = []
        for s in s_vals:
            if r >= s:
                continue
            st_t = data.get((graph, "ST", r, s))
            v3_t = data.get((graph, "V3", r, s))
            if st_t and st_t > 0 and v3_t and v3_t > 0:
                ss.append(s)
                speedups.append(st_t / v3_t)
            elif v3_t and v3_t > 0 and (st_t is None or st_t == -1):
                # ST failed, V3 OK → infinite speedup, cap at 1000
                ss.append(s)
                speedups.append(1000)
        if speedups:
            ax2.plot(ss, speedups, 'o-', label=f'r={r}', markersize=3)

    ax2.axhline(y=1, color='red', linestyle='--', alpha=0.5, label='V3=ST')
    ax2.set_xlabel("s")
    ax2.set_ylabel("Speedup (ST_time / V3_time)")
    ax2.set_yscale('log')
    ax2.legend(fontsize=8, ncol=3)
    ax2.grid(True, alpha=0.3)

    outf2 = f"bench_v3_{graph}_speedup.png"
    plt.savefig(outf2, dpi=150, bbox_inches='tight')
    print(f"  Saved {outf2}")
    plt.close()

    # --- V3 vs V3_NP speedup ---
    fig3, ax3 = plt.subplots(figsize=(12, 6))
    fig3.suptitle(f"{graph} — Private Cloud Speedup (V3 vs V3_NP)", fontsize=14)

    for r in r_vals:
        ss = []
        speedups = []
        for s in s_vals:
            if r >= s:
                continue
            np_t = data.get((graph, "V3_NP", r, s))
            v3_t = data.get((graph, "V3", r, s))
            if np_t and np_t > 0 and v3_t and v3_t > 0:
                ss.append(s)
                speedups.append(np_t / v3_t)
            elif v3_t and v3_t > 0 and (np_t is None or np_t == -1):
                ss.append(s)
                speedups.append(100)  # NP failed
        if speedups:
            ax3.plot(ss, speedups, 'o-', label=f'r={r}', markersize=3)

    ax3.axhline(y=1, color='red', linestyle='--', alpha=0.5)
    ax3.set_xlabel("s")
    ax3.set_ylabel("Speedup (V3_NP_time / V3_time)")
    ax3.set_yscale('log')
    ax3.legend(fontsize=8, ncol=3)
    ax3.grid(True, alpha=0.3)

    outf3 = f"bench_v3_{graph}_private_cloud.png"
    plt.savefig(outf3, dpi=150, bbox_inches='tight')
    print(f"  Saved {outf3}")
    plt.close()

print("\nDone.")
