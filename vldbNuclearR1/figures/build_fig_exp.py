"""
Figure: end-to-end performance vs. SOTA for R=1, 505 configurations on
three graphs (com-youtube, web-Stanford, web-it-2004).

Panel (a): per-graph wall-clock speedup curves over s (log-y).
Panel (b): cumulative distribution of speedups across all 505 configs.
Panel (c): peak-memory ratio (REF / Ours) per graph over s.

Source: benchmark_all_results.csv (Ours_ST vs REF_R1, status=OK).
"""
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['text.usetex']  = False
import csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict

OUT = Path(__file__).parent
CSV = OUT / "benchmark_all_results.csv"

rows = list(csv.DictReader(open(CSV)))
r1 = [r for r in rows if r['r'] == '1' and r['status'] == 'OK']
ours = {(r['graph'], int(r['s'])): r for r in r1 if r['algorithm'] == 'Ours_ST'}
ref  = {(r['graph'], int(r['s'])): r for r in r1 if r['algorithm'] == 'REF_R1'}
keys = sorted(set(ours) & set(ref))

# Per (graph, s): speedup and memory ratio
data = defaultdict(list)  # graph -> list of (s, speedup, mem_ratio)
all_speedup = []
for k in keys:
    g, s = k
    t_o = float(ours[k]['time_ms']);  t_r = float(ref[k]['time_ms'])
    m_o = float(ours[k]['memory_kB']); m_r = float(ref[k]['memory_kB'])
    if t_o <= 0 or m_o <= 0:
        continue
    sp = t_r / t_o
    mr = m_r / m_o
    data[g].append((s, sp, mr))
    all_speedup.append(sp)
for g in data:
    data[g].sort()

GRAPH_LABEL = {
    'com-youtube':  ('com-youtube',  '#1f77b4', 'o'),
    'web-Stanford': ('web-Stanford', '#2ca02c', 's'),
    'web-it-2004':  ('web-it-2004',  '#d62728', '^'),
}

fig = plt.figure(figsize=(15, 3.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.0], wspace=0.30)

# ============ Panel (a): speedup vs s ============
ax_a = fig.add_subplot(gs[0, 0])
for g in ['com-youtube', 'web-Stanford', 'web-it-2004']:
    xs = [d[0] for d in data[g]]
    ys = [d[1] for d in data[g]]
    label, col, mk = GRAPH_LABEL[g]
    ax_a.plot(xs, ys, color=col, marker=mk, markersize=3.5, linewidth=1.0,
              alpha=0.85, label=label)
ax_a.set_xscale('log')
ax_a.set_yscale('log')
ax_a.set_xlabel(r'parameter $s$', fontsize=10)
ax_a.set_ylabel(r'wall-clock speedup over SOTA', fontsize=10)
ax_a.set_title(r'(a) Speedup vs.\ $s$', fontsize=11)
ax_a.axhline(1.0, linestyle=':', color='#888', linewidth=0.8)
ax_a.grid(True, which='both', alpha=0.2)
ax_a.legend(fontsize=9, frameon=False, loc='upper left')

# ============ Panel (b): speedup CDF ============
ax_b = fig.add_subplot(gs[0, 1])
for g in ['com-youtube', 'web-Stanford', 'web-it-2004']:
    sps = sorted([d[1] for d in data[g]])
    ys = np.linspace(0, 1, len(sps))
    label, col, _ = GRAPH_LABEL[g]
    ax_b.plot(sps, ys, color=col, linewidth=1.6, label=label)
# overall
sps_all = sorted(all_speedup)
ax_b.plot(sps_all, np.linspace(0, 1, len(sps_all)),
          color='black', linewidth=1.0, linestyle='--',
          label=f'all 505 configs')
ax_b.set_xscale('log')
ax_b.set_xlabel(r'speedup over SOTA', fontsize=10)
ax_b.set_ylabel(r'cumulative fraction', fontsize=10)
ax_b.set_title(r'(b) Speedup distribution', fontsize=11)
ax_b.axvline(1.0, linestyle=':', color='#888', linewidth=0.8)
ax_b.grid(True, which='both', alpha=0.2)
ax_b.legend(fontsize=9, frameon=False, loc='lower right')

# Annotate geomean
import statistics
geo = statistics.geometric_mean(all_speedup)
ax_b.annotate(f'geomean $= {geo:.1f}\\times$', xy=(geo, 0.5),
              xytext=(geo*1.7, 0.18), fontsize=9,
              arrowprops=dict(arrowstyle='->', lw=0.8, color='#444'))

# ============ Panel (c): memory ratio ============
ax_c = fig.add_subplot(gs[0, 2])
for g in ['com-youtube', 'web-Stanford', 'web-it-2004']:
    xs = [d[0] for d in data[g]]
    ys = [d[2] for d in data[g]]
    label, col, mk = GRAPH_LABEL[g]
    ax_c.plot(xs, ys, color=col, marker=mk, markersize=3.5, linewidth=1.0,
              alpha=0.85, label=label)
ax_c.set_xscale('log')
ax_c.set_xlabel(r'parameter $s$', fontsize=10)
ax_c.set_ylabel(r'peak-RSS ratio (SOTA / Ours)', fontsize=10)
ax_c.set_title(r'(c) Memory reduction', fontsize=11)
ax_c.axhline(1.0, linestyle=':', color='#888', linewidth=0.8)
ax_c.grid(True, which='both', alpha=0.2)
ax_c.legend(fontsize=9, frameon=False, loc='upper right')

fig.savefig(OUT / 'fig_exp_endtoend.pdf', bbox_inches='tight')
fig.savefig(OUT / 'fig_exp_endtoend.png', dpi=170, bbox_inches='tight')
print(f"Saved fig_exp_endtoend to {OUT}")
print(f"  configs: {len(keys)}")
print(f"  speedup geomean = {geo:.2f}x, max = {max(all_speedup):.2f}x")
mem_ratios = [d[2] for g in data for d in data[g]]
print(f"  mem ratio max = {max(mem_ratios):.2f}x, geomean = {statistics.geometric_mean(mem_ratios):.2f}x")
