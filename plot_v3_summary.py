#!/usr/bin/env python3
"""Clean summary figure for V3 benchmark."""

import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

CSV = "/tmp/bench_v3_results.csv"
rows = list(csv.DictReader(open(CSV)))

# Parse
times = {}
ok_status = {}
for row in rows:
    g, r, s, algo, st = row["graph"], int(row["r"]), int(row["s"]), row["algo"], row["status"]
    wall = row.get("wall_ms", "")
    if st == "OK" and wall:
        times[(g, algo, r, s)] = float(wall)
    ok_status[(g, algo, r, s)] = st

graphs = ["dblp-core30", "com-dblp", "email-Eu-core", "com-youtube", "web-Stanford"]
graphs = [g for g in graphs if any(k[0] == g for k in times)]
colors = {"ST": "#e74c3c", "V3": "#27ae60", "V3_NP": "#3498db"}

# ========== Figure 1: Scalability frontier + OK counts ==========
fig, axes = plt.subplots(2, len(graphs), figsize=(4*len(graphs), 8),
                         gridspec_kw={'height_ratios': [2, 1]})
if len(graphs) == 1:
    axes = axes.reshape(2, 1)
fig.suptitle("V3 Private Cloud: Scalability Frontier", fontsize=15, fontweight='bold')

for idx, graph in enumerate(graphs):
    ax = axes[0, idx]
    for algo in ["ST", "V3", "V3_NP"]:
        s_vals = sorted(set(s for g, a, r, s in ok_status if g == graph and a == algo))
        frontier_s, frontier_r = [], []
        for sv in s_vals:
            max_r = max([r for r in range(3, sv)
                         if ok_status.get((graph, algo, r, sv)) == "OK"], default=0)
            if max_r > 0:
                frontier_s.append(sv)
                frontier_r.append(max_r)
        if frontier_s:
            ax.fill_between(frontier_s, frontier_r, alpha=0.12, color=colors[algo])
            ax.plot(frontier_s, frontier_r, '-', color=colors[algo],
                    linewidth=2, label=algo, alpha=0.9)

    ax.set_title(graph, fontsize=11, fontweight='bold')
    ax.set_xlabel("s", fontsize=10)
    if idx == 0:
        ax.set_ylabel("max r completed", fontsize=10)
        ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Bottom: OK count bars
    ax2 = axes[1, idx]
    algos = ["ST", "V3", "V3_NP"]
    counts = [sum(1 for k, v in ok_status.items() if k[0] == graph and k[1] == a and v == "OK")
              for a in algos]
    bars = ax2.bar(algos, counts, color=[colors[a] for a in algos], alpha=0.8)
    for bar, cnt in zip(bars, counts):
        if cnt > 0:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                     str(cnt), ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax2.set_ylabel("# OK" if idx == 0 else "")
    ax2.set_title("Solvable settings", fontsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("bench_v3_frontier.png", dpi=150, bbox_inches='tight')
print("Saved bench_v3_frontier.png")
plt.close()

# ========== Figure 2: Speedup analysis (dblp-core30 + com-dblp) ==========
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))
fig2.suptitle("V3 Speedup Analysis", fontsize=15, fontweight='bold')

# Panel 1: V3 vs ST speedup by r (dblp-core30)
ax = axes2[0]
for graph in ["dblp-core30", "com-dblp"]:
    all_r = sorted(set(r for g, a, r, s in times if g == graph and a == "V3"))
    r_speedups = []
    r_labels = []
    for rv in all_r:
        sps = []
        for s in range(rv+1, 200):
            st_t = times.get((graph, "ST", rv, s))
            v3_t = times.get((graph, "V3", rv, s))
            if st_t and v3_t and st_t > 0 and v3_t > 0:
                sps.append(st_t / v3_t)
        if sps:
            r_speedups.append(np.median(sps))
            r_labels.append(rv)
    if r_speedups:
        marker = 'o' if graph == "dblp-core30" else 's'
        ax.plot(r_labels, r_speedups, f'{marker}-', label=graph, linewidth=2, markersize=6)

ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel("r", fontsize=11)
ax.set_ylabel("Median speedup (ST/V3)", fontsize=11)
ax.set_title("V3 vs ST by r value", fontsize=12)
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 2: Private Cloud speedup by r (dblp-core30)
ax = axes2[1]
all_r = sorted(set(r for g, a, r, s in times if g == "dblp-core30" and a == "V3"))
r_labels, r_speedups = [], []
for rv in all_r:
    sps = []
    for s in range(rv+1, 200):
        np_t = times.get(("dblp-core30", "V3_NP", rv, s))
        v3_t = times.get(("dblp-core30", "V3", rv, s))
        if np_t and v3_t and np_t > 0 and v3_t > 0:
            sps.append(np_t / v3_t)
    if sps:
        r_labels.append(rv)
        r_speedups.append(np.median(sps))

ax.bar(range(len(r_labels)), r_speedups, color=colors["V3"], alpha=0.8)
ax.set_xticks(range(len(r_labels)))
ax.set_xticklabels(r_labels, fontsize=8)
ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel("r", fontsize=11)
ax.set_ylabel("Median speedup (V3_NP / V3)", fontsize=11)
ax.set_title("Private Cloud effect (dblp-core30)", fontsize=12)
ax.set_yscale('log')
ax.grid(True, alpha=0.3, axis='y')

# Panel 3: Wall time heatmap for V3 on dblp-core30
ax = axes2[2]
g = "dblp-core30"
all_r = sorted(set(r for gg, a, r, s in times if gg == g and a == "V3"))
all_s = sorted(set(s for gg, a, r, s in times if gg == g and a == "V3"))
mat = np.full((len(all_r), len(all_s)), np.nan)
for i, r in enumerate(all_r):
    for j, s in enumerate(all_s):
        t = times.get((g, "V3", r, s))
        if t and t > 0:
            mat[i, j] = np.log10(t)

im = ax.imshow(mat, aspect='auto', origin='lower', cmap='YlOrRd', vmin=2.5, vmax=6)
r_step = max(1, len(all_r) // 8)
s_step = max(1, len(all_s) // 8)
ax.set_yticks(range(0, len(all_r), r_step))
ax.set_yticklabels([all_r[i] for i in range(0, len(all_r), r_step)], fontsize=8)
ax.set_xticks(range(0, len(all_s), s_step))
ax.set_xticklabels([all_s[i] for i in range(0, len(all_s), s_step)], fontsize=8)
ax.set_xlabel("s", fontsize=11)
ax.set_ylabel("r", fontsize=11)
ax.set_title("V3 wall time (dblp-core30)", fontsize=12)
cbar = fig2.colorbar(im, ax=ax, shrink=0.8)
cbar.set_ticks([3, 4, 5, 6])
cbar.set_ticklabels(['1s', '10s', '100s', '1000s'])

plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig("bench_v3_speedup.png", dpi=150, bbox_inches='tight')
print("Saved bench_v3_speedup.png")
plt.close()

print("Done.")
