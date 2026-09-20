#!/usr/bin/env python3
"""Generate the paper's figures (figures/*.pdf) from the evidence JSON files of research/r1_skyline_index_20260918.
Run from anywhere; each figure function skips silently when its evidence is not there yet."""
import json
import statistics
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
EV = HERE.parent / 'research' / 'r1_skyline_index_20260918'
OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 7.5, 'axes.labelsize': 7.5, 'legend.fontsize': 6.5, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                     'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
NAMES = {'soc-pokec-relationships': 'soc-pokec', 'com-amazon.ungraph': 'com-amazon'}
COLORS = {'laptop': '#1f5f8b', 'tods1': '#b5541c', 'tods2': '#3b7a3b'}

def load(name):
    p = EV / name
    return json.loads(p.read_text()) if p.exists() else None

def merged_rows():
    best = {}
    for where, f in [('tods1', 'tods1.json'), ('tods2', 'tods2.json'), ('laptop', 'final.json'), ('laptop', 'more.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' not in r: continue
            key = (r['result']['n'], r['result']['m'])
            if key not in best: best[key] = (where, NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem), r['result'])
    return list(best.values())

def fig_regimes():
    """Community queries at the three levels: listing time against community size (left), locating time (right)."""
    rows = merged_rows()
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.1), gridspec_kw={'width_ratios': [1.35, 1]})
    marks = {'own': 'o', 'half': 's', 'root': '^'}
    for where in ('laptop', 'tods1', 'tods2'):
        for reg, mk in marks.items():
            xs = [x[f'{reg}_output'] for w, g, x in rows if w == where]; ys = [x[f'explicit_{reg}_ns'] for w, g, x in rows if w == where]
            a.scatter(xs, ys, s=11, marker=mk, facecolors='none' if reg != 'own' else COLORS[where], edgecolors=COLORS[where], linewidths=0.7)
    lo, hi = 1, 3e6
    for slope, lab in ((0.05, '0.05 ns per vertex'), (0.5, '0.5 ns per vertex')):
        a.plot([lo, hi], [slope * lo + 8, slope * hi + 8], color='0.6', lw=0.6, ls='--'); a.text(hi, slope * hi + 8, ' ' + lab, fontsize=6, color='0.4', va='center')
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('community size (vertices, mean over the queries)'); a.set_ylabel('listing time (ns)')
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker=mk, color='0.3', ls='', markerfacecolor='0.3' if reg == 'own' else 'none', markersize=4, label={'own': 'own level', 'half': 'half level', 'root': 'k = 1'}[reg]) for reg, mk in marks.items()]
    handles += [Line2D([], [], marker='o', color=c, ls='', markersize=4, label=w) for w, c in COLORS.items()]
    a.legend(handles=handles, ncol=2, loc='upper left', handletextpad=0.2, columnspacing=0.8)
    a.set_xlim(lo, hi * 8)
    data = [[x[f'ptr_{reg}_ns'] for w, g, x in rows] for reg in marks]
    b.boxplot(data, tick_labels=['own level', 'half level', 'k = 1'], widths=0.5, medianprops={'color': '#b5541c'}, flierprops={'markersize': 3})
    for i, d in enumerate(data):
        b.scatter([i + 1 + (j % 5 - 2) * 0.05 for j in range(len(d))], d, s=6, color='0.35', zorder=3)
    b.set_ylabel('locating time (ns)'); b.set_ylim(0, None)
    fig.tight_layout(w_pad=1.5); fig.savefig(OUT / 'fig_regimes.pdf'); plt.close(fig)

def profile_runs(tag):
    rec = load(f'{tag}.json')
    if rec is None: return []
    return [(NAMES.get(Path(r['index']).stem, Path(r['index']).stem), r['result']) for r in rec['runs'] if 'result' in r]

def fig_profile(tag='profile_tods2'):
    """Own-level latency by clique size (stratified workload) and by community size (deciles of the fixed workload)."""
    runs = [(g, x) for g, x in profile_runs(tag) if '_p' not in g]
    if not runs: return
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.0))
    cmap = plt.get_cmap('tab10')
    for i, (g, x) in enumerate(runs):
        bs = x['by_size']
        a.plot([e['s'] for e in bs], [e['locate_ns'] for e in bs], marker='o', ms=2.2, lw=0.8, color=cmap(i), label=g)
        b.plot([e['s'] for e in bs], [e['list_ns'] / max(e['output'], 1) for e in bs], marker='o', ms=2.2, lw=0.8, color=cmap(i))
        dc = x['own_by_output']
        c.plot([e['output'] for e in dc], [e['list_ns'] for e in dc], marker='o', ms=2.2, lw=0.8, color=cmap(i))
    a.set_xscale('log'); a.set_xlabel('clique size s'); a.set_ylabel('locating time (ns)'); a.set_ylim(0, None)
    a.legend(ncol=1, loc='upper right', handlelength=1.2, handletextpad=0.3, labelspacing=0.2)
    b.set_xscale('log'); b.set_yscale('log'); b.set_xlabel('clique size s'); b.set_ylabel('listing time per vertex (ns)')
    c.set_xscale('log'); c.set_yscale('log'); c.set_xlabel('community size (vertices, decile mean)'); c.set_ylabel('listing time (ns)')
    fig.tight_layout(w_pad=1.2); fig.savefig(OUT / 'fig_profile.pdf'); plt.close(fig)

def fig_scale(tag='scale_tods2', full='tods2.json'):
    """Vertex-induced samples at 20 to 100 percent: index bytes, S trees bytes, build time, own-level latencies."""
    rec = load(f'{tag}.json'); fullrec = load(full)
    if rec is None or fullrec is None: return
    series = {}
    for r in rec['runs']:
        if 'result' not in r: continue
        stem = Path(r['graph']).stem; g, p = stem.rsplit('_p', 1); series.setdefault(g, {})[int(p)] = (r['result'], r.get('peak_rss_bytes'))
    for r in fullrec['runs']:
        if 'result' in r and Path(r['graph']).stem in series: series[Path(r['graph']).stem][100] = (r['result'], r.get('peak_rss_bytes'))
    fig, axes = plt.subplots(1, 4, figsize=(7.0, 1.9))
    cmap = plt.get_cmap('tab10')
    for i, (g, pts) in enumerate(sorted(series.items())):
        ps = sorted(pts); xs = [pts[p][0]['n'] / 1e6 for p in ps]
        axes[0].plot(xs, [pts[p][0]['bytes_total'] / 1048576 for p in ps], marker='o', ms=2.5, lw=0.8, color=cmap(i), label=f'{g}, index')
        axes[0].plot(xs, [pts[p][0]['baseline_vertex_bytes'] / 1048576 for p in ps], marker='s', ms=2.5, lw=0.8, ls='--', color=cmap(i), label=f'{g}, S trees')
        axes[1].plot(xs, [(pts[p][0]['ti_ms'] + pts[p][0]['build_ms'] + pts[p][0]['compact_ms']) / 1000 for p in ps], marker='o', ms=2.5, lw=0.8, color=cmap(i), label=g)
        axes[2].plot(xs, [pts[p][0]['ptr_own_ns'] for p in ps], marker='o', ms=2.5, lw=0.8, color=cmap(i), label=g)
        axes[3].plot(xs, [pts[p][0]['explicit_own_ns'] / pts[p][0]['own_output'] for p in ps], marker='o', ms=2.5, lw=0.8, color=cmap(i), label=g)
    for ax, lab in zip(axes, ('bytes (MB)', 'build time (s)', 'locating time (ns)', 'listing time per vertex (ns)')):
        ax.set_xlabel('vertices (millions)'); ax.set_ylabel(lab); ax.set_ylim(0, None); ax.set_xlim(0, None)
    axes[0].set_yscale('log'); axes[0].set_ylim(None, None); axes[0].legend(handlelength=1.5, labelspacing=0.2)
    fig.tight_layout(w_pad=1.0); fig.savefig(OUT / 'fig_scale.pdf'); plt.close(fig)

if __name__ == '__main__':
    fig_regimes(); fig_profile(); fig_scale()
    print('figures written to', OUT)
