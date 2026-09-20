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

def profile_runs(tag):
    rec = load(f'{tag}.json')
    if rec is None: return []
    return [(NAMES.get(Path(r['index']).stem, Path(r['index']).stem), r['result']) for r in rec['runs'] if 'result' in r]

def all_profiles(samples=False):
    """(machine, graph, result) over the profile records; samples excluded unless asked; one row per graph, servers first."""
    seen = {}
    for where in ('tods2', 'tods1', 'laptop'):
        for g, x in profile_runs(f'profile_{where}'):
            if ('_p' in g) != samples: continue
            seen.setdefault((x['n'], x['s_max']), (where, g, x))
    return list(seen.values())

def fig_regimes():
    """Community queries at the three levels, index against S trees: listing time against answer size (left), locating time (right)."""
    rows = all_profiles()
    if not rows: return
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.1), gridspec_kw={'width_ratios': [1.35, 1]})
    marks = {'own': 'o', 'half': 's', 'root': '^'}; labels = {'own': 'own level', 'half': 'half level', 'root': 'k = 1'}
    for reg, mk in marks.items():
        xs = [x['regimes'][reg]['output'] for w, g, x in rows]
        a.scatter(xs, [x['regimes'][reg]['list_ns'] for w, g, x in rows], s=12, marker=mk, color='#1f5f8b', linewidths=0.6, zorder=3)
        a.scatter(xs, [x['regimes'][reg]['st_list_ns'] for w, g, x in rows], s=14, marker=mk, facecolors='none', edgecolors='#b5541c', linewidths=0.7, zorder=2)
    lo, hi = 10, 4e6
    for slope, lab in ((0.05, '0.05 ns per vertex'), (0.5, '0.5 ns per vertex')):
        a.plot([lo, hi], [slope * lo + 8, slope * hi + 8], color='0.6', lw=0.6, ls='--'); a.text(hi, slope * hi + 8, ' ' + lab, fontsize=6, color='0.4', va='center')
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel('answer size (vertices, mean over the queries)'); a.set_ylabel('listing time (ns)'); a.set_xlim(lo, hi * 8)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker=mk, color='0.3', ls='', markersize=4, label=labels[reg]) for reg, mk in marks.items()]
    handles += [Line2D([], [], marker='o', color='#1f5f8b', ls='', markersize=4, label='chain index'), Line2D([], [], marker='o', color='#b5541c', markerfacecolor='none', ls='', markersize=4, label='S trees')]
    a.legend(handles=handles, ncol=2, loc='upper left', handletextpad=0.2, columnspacing=0.8)
    data = []; ticks = []
    for reg in marks:
        data.append([x['regimes'][reg]['locate_ns'] for w, g, x in rows]); data.append([x['regimes'][reg]['st_locate_ns'] for w, g, x in rows]); ticks.append(labels[reg])
    pos = [1, 1.7, 3, 3.7, 5, 5.7]
    bp = b.boxplot(data, positions=pos, widths=0.55, medianprops={'color': '0.2'}, flierprops={'markersize': 2.5}, patch_artist=True)
    for i, box in enumerate(bp['boxes']): box.set(facecolor='#1f5f8b' if i % 2 == 0 else 'white', edgecolor='#1f5f8b' if i % 2 == 0 else '#b5541c', alpha=0.85)
    b.set_xticks([1.35, 3.35, 5.35]); b.set_xticklabels(ticks); b.set_ylabel('locating time (ns)'); b.set_yscale('log')
    b.legend(handles=[Line2D([], [], color='#1f5f8b', lw=4, label='chain index'), Line2D([], [], color='#b5541c', lw=1.2, label='S trees')], loc='upper right')
    fig.tight_layout(w_pad=1.5); fig.savefig(OUT / 'fig_regimes.pdf'); plt.close(fig)

def fig_profile():
    """Own-level latency by clique size (stratified workload) and by answer size (deciles), index (solid) against S trees (dashed)."""
    rows = [(w, g, x) for w, g, x in all_profiles() if w != 'laptop'] or all_profiles()
    if not rows: return
    rows = sorted(rows, key=lambda t: -t[2]['n'])[:4]
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.0))
    cmap = plt.get_cmap('tab10')
    for i, (w, g, x) in enumerate(rows):
        bs = x['by_size']; col = cmap(i)
        a.plot([e['s'] for e in bs], [e['locate_ns'] for e in bs], marker='o', ms=2, lw=0.8, color=col, label=g)
        a.plot([e['s'] for e in bs], [e['st_locate_ns'] for e in bs], marker='o', ms=2, lw=0.8, ls='--', color=col, markerfacecolor='none')
        b.plot([e['s'] for e in bs], [e['list_ns'] for e in bs], marker='o', ms=2, lw=0.8, color=col)
        b.plot([e['s'] for e in bs], [e['st_list_ns'] for e in bs], marker='o', ms=2, lw=0.8, ls='--', color=col, markerfacecolor='none')
        dc = x['own_by_output']
        c.plot([e['output'] for e in dc], [e['list_ns'] for e in dc], marker='o', ms=2, lw=0.8, color=col)
        c.plot([e['output'] for e in dc], [e['st_list_ns'] for e in dc], marker='o', ms=2, lw=0.8, ls='--', color=col, markerfacecolor='none')
    from matplotlib.lines import Line2D
    a.set_xscale('log'); a.set_xlabel('clique size s'); a.set_ylabel('locating time (ns)'); a.set_ylim(0, None)
    a.legend(ncol=1, loc='upper left', handlelength=1.2, handletextpad=0.3, labelspacing=0.15, fontsize=5.5)
    b.set_xscale('log'); b.set_yscale('log'); b.set_xlabel('clique size s'); b.set_ylabel('listing time (ns)')
    b.legend(handles=[Line2D([], [], color='0.3', lw=0.9, label='chain index'), Line2D([], [], color='0.3', lw=0.9, ls='--', label='S trees')], loc='upper right')
    c.set_xscale('log'); c.set_yscale('log'); c.set_xlabel('answer size (vertices, decile mean)'); c.set_ylabel('listing time (ns)')
    fig.tight_layout(w_pad=1.2); fig.savefig(OUT / 'fig_profile.pdf'); plt.close(fig)

def prior_total(where, g):
    """CND summed build and peel time (s) over all sizes for (machine, graph), from the prior records."""
    d = {'laptop': 'prior', 'tods2': 'prior/tods2', 'tods1': 'prior/tods1'}[where]
    p = EV / d / f'prior_original_{g}.json'
    if not p.exists(): p = EV / d / f"prior_original_{ {v: k for k, v in NAMES.items()}.get(g, g) }.json"
    if not p.exists(): return None
    og = json.loads(p.read_text()); return og['total_inproc_ms'] / 1000 if 'total_inproc_ms' in og else None

def fig_scale(tag='scale_tods2', full='tods2.json'):
    """Vertex-induced samples at 20 to 100 percent: bytes (index, S trees), build time (index, CND all sizes), own-level latencies (index, S trees)."""
    rec = load(f'{tag}.json'); fullrec = load(full)
    if rec is None or fullrec is None: return
    series = {}
    for r in rec['runs']:
        if 'result' not in r: continue
        stem = Path(r['graph']).stem; g, p = stem.rsplit('_p', 1); series.setdefault(g, {})[int(p)] = r['result']
    for r in fullrec['runs']:
        if 'result' in r and Path(r['graph']).stem in series: series[Path(r['graph']).stem][100] = r['result']
    prof = {g: x for w, g, x in all_profiles(samples=True)}; prof.update({g: x for w, g, x in all_profiles() if w == 'tods2'})
    fig, axes = plt.subplots(1, 4, figsize=(7.0, 1.9))
    cmap = plt.get_cmap('tab10')
    from matplotlib.lines import Line2D
    for i, (g, pts) in enumerate(sorted(series.items())):
        ps = sorted(pts); xs = [pts[p]['n'] / 1e6 for p in ps]; col = cmap(i)
        axes[0].plot(xs, [pts[p]['bytes_total'] / 1048576 for p in ps], marker='o', ms=2.5, lw=0.8, color=col, label=g)
        axes[0].plot(xs, [pts[p]['baseline_vertex_bytes'] / 1048576 for p in ps], marker='o', ms=2.5, lw=0.8, ls='--', color=col, markerfacecolor='none')
        axes[1].plot(xs, [(pts[p]['ti_ms'] + pts[p]['build_ms'] + pts[p]['compact_ms']) / 1000 for p in ps], marker='o', ms=2.5, lw=0.8, color=col)
        cnd = [(pts[p]['n'] / 1e6, prior_total('tods2', g if p == 100 else f'{g}_p{p}')) for p in ps]; cnd = [(x, y) for x, y in cnd if y is not None]
        if cnd: axes[1].plot([x for x, y in cnd], [y for x, y in cnd], marker='o', ms=2.5, lw=0.8, ls='--', color=col, markerfacecolor='none')
        names = [g if p == 100 else f'{g}_p{p}' for p in ps]
        have = [(pts[p]['n'] / 1e6, prof[nm]) for p, nm in zip(ps, names) if nm in prof]
        if have:
            axes[2].plot([x for x, y in have], [y['regimes']['own']['locate_ns'] for x, y in have], marker='o', ms=2.5, lw=0.8, color=col)
            axes[2].plot([x for x, y in have], [y['regimes']['own']['st_locate_ns'] for x, y in have], marker='o', ms=2.5, lw=0.8, ls='--', color=col, markerfacecolor='none')
            axes[3].plot([x for x, y in have], [y['regimes']['own']['list_ns'] / y['regimes']['own']['output'] for x, y in have], marker='o', ms=2.5, lw=0.8, color=col)
            axes[3].plot([x for x, y in have], [y['regimes']['own']['st_list_ns'] / y['regimes']['own']['output'] for x, y in have], marker='o', ms=2.5, lw=0.8, ls='--', color=col, markerfacecolor='none')
    for ax, lab in zip(axes, ('bytes (MB)', 'build time, all sizes (s)', 'locating time (ns)', 'listing time per vertex (ns)')):
        ax.set_xlabel('vertices (millions)'); ax.set_ylabel(lab); ax.set_xlim(0, None)
    axes[0].set_yscale('log'); axes[1].set_yscale('log'); axes[2].set_ylim(0, None); axes[3].set_ylim(0, None)
    axes[0].legend(handlelength=1.5, labelspacing=0.2, loc='lower right')
    axes[1].legend(handles=[Line2D([], [], color='0.3', lw=0.9, label='chain index'), Line2D([], [], color='0.3', lw=0.9, ls='--', label='S trees / CND')], loc='lower right')
    fig.tight_layout(w_pad=1.0); fig.savefig(OUT / 'fig_scale.pdf'); plt.close(fig)

if __name__ == '__main__':
    fig_regimes(); fig_profile(); fig_scale()
    print('figures written to', OUT)
