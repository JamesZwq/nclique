#!/usr/bin/env python3
"""Plot benchmark results for all graphs, all r values.
Usage: python3 plot_bench_all.py [csv_file]
  Default csv: benchmark_combined.csv
  Outputs: bench_<graph>_all.png per graph
"""
import csv
import sys
import os
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 10

TIMEOUT_SEC = 3600
csvfile = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), 'benchmark_combined.csv')
outdir = os.path.dirname(os.path.abspath(csvfile))

data = {}
with open(csvfile) as f:
    reader = csv.DictReader(f)
    for row in reader:
        g = row['graph']
        r = int(row['r'])
        s = int(row['s'])
        algo = row['algorithm']
        st = row['status']
        if st == 'OK':
            t = float(row['time_ms']) / 1000
            mem = float(row['memory_kB']) / 1024 if row['memory_kB'] not in ('N/A', '') else None
        elif st == 'TIMEOUT':
            t = TIMEOUT_SEC
            mem = None
        else:
            continue
        key = (g, r, algo)
        data.setdefault(key, {'s': [], 't': [], 'm': []})
        if s not in data[key]['s']:
            data[key]['s'].append(s)
            data[key]['t'].append(t)
            data[key]['m'].append(mem)

for k in data:
    pairs = sorted(zip(data[k]['s'], data[k]['t'], data[k]['m']))
    data[k]['s'] = [p[0] for p in pairs]
    data[k]['t'] = [p[1] for p in pairs]
    data[k]['m'] = [p[2] for p in pairs]

algo_config = {
    1: [('Ours_ST', 'Ours', '#2196F3', 'o'),
        ('REF_R1', 'REF', '#FF5722', 's')],
    2: [('Ours_DCLP', 'Ours', '#2196F3', 'o'),
        ('REF_R2', 'REF', '#FF5722', 's')],
    3: [('Ours_V20', 'Ours', '#2196F3', 'o'),
        ('REF_R3', 'REF', '#FF5722', 's')],
    4: [('Ours_V20', 'Ours', '#2196F3', 'o'),
        ('REF_R3', 'REF', '#FF5722', 's')],
    5: [('Ours_V20', 'Ours', '#2196F3', 'o'),
        ('REF_R3', 'REF', '#FF5722', 's')],
}

graphs_info = {
    'com-dblp': 'com-dblp (317K V, 1.05M E)',
    'com-youtube': 'com-youtube (1.13M V, 2.99M E)',
    'web-Stanford': 'web-Stanford (282K V, 2.31M E)',
    'web-it-2004': 'web-it-2004 (509K V, 7.18M E)',
}

graphs_found = sorted(set(k[0] for k in data))

for gname in graphs_found:
    gtitle = graphs_info.get(gname, gname)
    r_values = sorted(set(k[1] for k in data if k[0] == gname))
    if not r_values:
        continue
    ncols = len(r_values)
    fig, axes = plt.subplots(2, ncols, figsize=(5 * ncols, 8))
    if ncols == 1:
        axes = axes.reshape(2, 1)

    for ri, r in enumerate(r_values):
        algos = algo_config.get(r, algo_config[3])
        ax_t = axes[0][ri]
        ax_m = axes[1][ri]

        for akey, alabel, color, marker in algos:
            d = data.get((gname, r, akey))
            if not d or len(d['s']) < 2:
                continue
            ax_t.plot(d['s'], d['t'], f'{marker}-', color=color, label=alabel,
                      linewidth=1.5, markersize=4, alpha=0.8)
            mems = [m for m in d['m'] if m is not None]
            ss_m = [s for s, m in zip(d['s'], d['m']) if m is not None]
            if mems:
                ax_m.plot(ss_m, mems, f'{marker}-', color=color, label=alabel,
                          linewidth=1.5, markersize=4, alpha=0.8)

        ax_t.axhline(y=TIMEOUT_SEC, color='gray', linestyle='--', alpha=0.5, label='Timeout (1h)')
        ax_t.set_xlabel('s')
        if ri == 0:
            ax_t.set_ylabel('Time (seconds)')
        ax_t.set_title(f'r = {r}')
        ax_t.legend(fontsize=8, loc='best')
        ax_t.grid(True, alpha=0.3)
        ax_t.set_yscale('log')

        ax_m.set_xlabel('s')
        if ri == 0:
            ax_m.set_ylabel('Memory (MB)')
        ax_m.set_title(f'r = {r}')
        ax_m.legend(fontsize=8, loc='best')
        ax_m.grid(True, alpha=0.3)

    fig.suptitle(gtitle, fontsize=14, fontweight='bold')
    plt.tight_layout()
    fname = os.path.join(outdir, f'bench_{gname}_all.png')
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {fname}')
