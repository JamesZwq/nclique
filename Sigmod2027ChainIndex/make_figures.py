#!/usr/bin/env python3
"""Generate the paper's experiment figures (figures/*.pdf) from the evidence JSON files of
research/r1_skyline_index_20260918.  Run from anywhere; each figure skips silently when its evidence is missing.

House style: monochrome; ChainIndex solid black with filled markers, the baseline (STrees or CND) dotted grey with
hollow markers; every figure is drawn at exactly the width it is printed at (a figure* is \\textwidth = 506.3pt),
so the 8pt type below is the type size on the page; serif type matching the paper (Linux Libertine)."""
import json
import statistics
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# the paper's font, installed as user fonts on this machine; registered explicitly since matplotlib's cache skips them
from matplotlib import font_manager as _fm
import glob as _glob
for _f in _glob.glob(str(Path.home() / 'Library' / 'Fonts' / 'LinLibertine_R*.ttf')):
    try: _fm.fontManager.addfont(_f)
    except Exception: pass

HERE = Path(__file__).resolve().parent
EV = HERE.parent / 'research' / 'r1_skyline_index_20260918'
OUT = HERE / 'figures'
OUT.mkdir(exist_ok=True)
PT = 1 / 72.0
TEXT_WIDTH = 506.295 * PT
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Linux Libertine', 'Linux Libertine O', 'Times New Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'custom', 'mathtext.rm': 'Linux Libertine', 'mathtext.it': 'Linux Libertine:italic',
    'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7.4, 'ytick.labelsize': 7.4, 'legend.fontsize': 7.2,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
    'xtick.major.size': 2.2, 'ytick.major.size': 2.2, 'xtick.minor.size': 1.2, 'ytick.minor.size': 1.2,
    'legend.frameon': False, 'legend.handlelength': 1.6, 'legend.handletextpad': 0.4, 'legend.labelspacing': 0.25,
    'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42,
})
NAMES = {'soc-pokec-relationships': 'soc-pokec', 'com-amazon.ungraph': 'com-amazon'}
OURS, BASE = 'black', '0.5'
OURS_KW = dict(color=OURS, lw=0.9, ls='-')
BASE_KW = dict(color=BASE, lw=0.9, ls=(0, (1.2, 1.4)), markerfacecolor='white')
MARKS = {'web-BerkStan': 'o', 'cit-Patents': 's'}
REP = ['web-BerkStan', 'cit-Patents']                      # the representative pair of Exp-5 and Exp-7

def load(name):
    p = EV / name
    return json.loads(p.read_text()) if p.exists() else None

def profile_runs(tag):
    rec = load(f'{tag}.json')
    if rec is None: return []
    return [(NAMES.get(Path(r['index']).stem, Path(r['index']).stem), r['result']) for r in rec['runs'] if 'result' in r]

def all_profiles(samples=False, pairs=False):
    """(machine, graph, result) over the profile records; samples excluded unless asked; one row per graph, servers
    first, or every (graph, machine) pair when pairs=True (the set the text's level medians are computed on)."""
    seen = {}
    for where in ('tods2', 'tods1', 'laptop'):
        for g, x in profile_runs(f'profile_{where}'):
            if ('_p' in g) != samples: continue
            seen.setdefault((x['n'], x['s_max']) if not pairs else (where, g), (where, g, x))
    return list(seen.values())

def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf'); plt.close(fig); print('  wrote figures/%s.pdf' % name)

def legend_pair(ax, base_label='STrees', loc='best', **kw):
    ax.legend(handles=[Line2D([], [], marker='o', ms=3.2, **OURS_KW, label='ChainIndex'),
                       Line2D([], [], marker='o', ms=3.2, **BASE_KW, label=base_label)], loc=loc, **kw)

# ------------------------------------------------------------------------------- Exp-4: by level ----
def fig_regimes():
    """Left: listing time against the mean answer size, one point per (graph, level), both indexes.
    Right: locating time by level, one point per graph, both indexes."""
    rows = all_profiles(pairs=True)
    if not rows: return
    # 2026-09-22: three panels; (b) replaces the per-graph query table: the listing-time ratio against vertices per range
    fig, (a, c, b) = plt.subplots(1, 3, figsize=(TEXT_WIDTH, 1.55), gridspec_kw={'width_ratios': [1.3, 1, 1]})
    marks = {'own': 'o', 'half': 's', 'root': '^'}; labels = {'own': 'own level', 'half': 'half level', 'root': 'k = 1'}
    for reg, mk in marks.items():
        xs = [x['regimes'][reg]['output'] for w, g, x in rows]
        a.scatter(xs, [x['regimes'][reg]['list_ns'] for w, g, x in rows], s=11, marker=mk, color=OURS, linewidths=0.5, zorder=3)
        a.scatter(xs, [x['regimes'][reg]['st_list_ns'] for w, g, x in rows], s=13, marker=mk, facecolors='white', edgecolors=BASE, linewidths=0.7, zorder=2)
    lo, hi = 10, 4e6
    for slope, lab in ((0.05, '0.05 ns per vertex'), (0.5, '0.5 ns per vertex')):
        a.plot([lo, hi], [slope * lo + 8, slope * hi + 8], color='0.7', lw=0.6, ls=(0, (2, 1.5)))
        a.text(hi * 1.15, slope * hi + 8, lab, fontsize=6.6, color='0.35', va='center')
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(lo, hi * 12); a.set_ylim(5, 5e6)
    a.set_xlabel('mean answer size (vertices)'); a.set_ylabel('listing time (ns)')
    handles = [Line2D([], [], marker=mk, color='0.25', ls='', markersize=3.6, label=labels[reg]) for reg, mk in marks.items()]
    handles += [Line2D([], [], marker='o', color=OURS, ls='', markersize=3.6, label='ChainIndex'),
                Line2D([], [], marker='o', color=BASE, markerfacecolor='white', ls='', markersize=3.6, label='STrees')]
    a.legend(handles=handles, ncol=2, loc='upper left', columnspacing=0.9, handletextpad=0.2)
    # middle: the listing time of STrees divided by that of ChainIndex against the mean vertices per range
    for reg, mk in marks.items():
        pts = [(x['regimes'][reg]['output'] / max(x['regimes'][reg]['ranges'], 1), x['regimes'][reg]['st_list_ns'] / x['regimes'][reg]['list_ns']) for w, g, x in rows]
        c.scatter([p[0] for p in pts], [p[1] for p in pts], s=11, marker=mk, color=OURS, linewidths=0.5, zorder=3)
    c.axhline(1, color='0.7', lw=0.6, ls=(0, (2, 1.5)))
    c.text(1.2e4, 0.17, 'ChainIndex faster above 1', fontsize=6.2, color='0.35', ha='right', va='bottom')
    c.set_xscale('log'); c.set_yscale('log'); c.set_xlabel('mean vertices per range'); c.set_ylabel('listing time ratio (STrees / ChainIndex)')
    c.set_xlim(1, 2e4); c.set_ylim(0.15, 8)
    c.set_yticks([0.2, 0.5, 1, 2, 5]); c.set_yticklabels(['0.2', '0.5', '1', '2', '5'])
    # right: strip plot of the locating times, one point per graph, the two indexes side by side
    import numpy as np
    rng = np.random.default_rng(3)
    for i, reg in enumerate(marks):
        ours = [x['regimes'][reg]['locate_ns'] for w, g, x in rows]; st = [x['regimes'][reg]['st_locate_ns'] for w, g, x in rows]
        jit = rng.uniform(-0.09, 0.09, len(ours))
        b.scatter(i - 0.16 + jit, ours, s=9, marker='o', color=OURS, linewidths=0.5, zorder=3)
        b.scatter(i + 0.16 + jit, st, s=11, marker='o', facecolors='white', edgecolors=BASE, linewidths=0.7, zorder=2)
        top = max(max(ours), max(st))
        for off, vals, col in ((-0.16, ours, OURS), (0.16, st, BASE)):
            med = statistics.median(vals); b.plot([i + off - 0.15, i + off + 0.15], [med, med], color=col, lw=1.0, zorder=4)
            b.text(i + off, top * 1.7, f'{med:.0f}', ha='center', va='bottom', fontsize=6.8, color=col, zorder=5)
    b.text(2.55, 4000 * 0.55, 'medians', ha='right', va='top', fontsize=6.6, color='0.35')
    b.set_xticks(range(3)); b.set_xticklabels([labels[r] for r in marks]); b.set_xlim(-0.6, 2.6)
    b.set_yscale('log'); b.set_ylabel('locating time (ns)'); b.set_ylim(1.5, 6000)
    legend_pair(b, loc='upper left')
    for ax, t in zip((a, c, b), ('(a)', '(b)', '(c)')): ax.set_title(t, loc='left', fontsize=8, pad=3)
    fig.subplots_adjust(left=0.07, right=0.995, bottom=0.23, top=0.89, wspace=0.4)
    save(fig, 'fig_regimes')

# ------------------------------------------------------------- Exp-8: CND once per size against one build ----
COL_WIDTH = 240.96 * PT

def fig_prior():
    """One point per (graph, machine) pair: the summed build and peel time of CND over all sizes against the build
    time of ChainIndex on the same machine, log-log, with the 10x/100x/1000x lines.  Replaces the per-pair table."""
    ours = {}
    for where, f in [('laptop', 'final.json'), ('laptop', 'more.json'), ('tods1', 'tods1.json'), ('tods2', 'tods2.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' in r: ours.setdefault((where, NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem)), r['result'])
    pts = []
    for where, d in [('laptop', 'prior'), ('tods2', 'prior/tods2'), ('tods1', 'prior/tods1')]:
        for p in sorted((EV / d).glob('prior_original_*.json')):
            og = json.loads(p.read_text()); g = NAMES.get(p.stem[len('prior_original_'):], p.stem[len('prior_original_'):])
            if 'total_inproc_ms' in og and (where, g) in ours:
                o = ours[(where, g)]; pts.append((where, g, (o['ti_ms'] + o['build_ms'] + o['compact_ms']) / 1000, og['total_inproc_ms'] / 1000, og['sizes_ok']))
    if not pts: return
    fig, ax = plt.subplots(figsize=(COL_WIDTH, 1.3))
    lo, hi = 5e-3, 150
    for f, lab in ((1, '1x'), (10, '10x'), (100, '100x'), (1000, '1000x')):
        ax.plot([lo, hi], [lo * f, hi * f], color='0.75', lw=0.6, ls=(0, (2, 1.5)), zorder=1)
        ax.text(hi * 1.15, hi * f, lab, fontsize=6.4, color='0.4', va='center')
    mk = {'laptop': 'o', 'tods1': 's', 'tods2': '^'}; lab = {'laptop': 'laptop', 'tods1': 'server 1', 'tods2': 'server 2'}
    for where in mk:
        sel = [p for p in pts if p[0] == where]
        if sel: ax.scatter([p[2] for p in sel], [p[3] for p in sel], s=12, marker=mk[where], color=OURS, linewidths=0.5, zorder=3, label=lab[where])
    done = set()
    for where, g, x, y, k in sorted(pts, key=lambda p: p[0] != 'tods1'):
        if g in ('web-uk-2005', 'com-amazon', 'web-BerkStan') and g not in done:
            done.add(g); ax.annotate(f'{g}, {k} sizes', (x, y), textcoords='offset points', xytext=(-5, 7), ha='right', fontsize=6.2, color='0.25')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(lo, 1500); ax.set_ylim(0.05, 4e5)
    ax.set_xlabel('ChainIndex, one build for every size (s)'); ax.set_ylabel('CND, one run per size (s)')
    ax.legend(loc='lower right', handletextpad=0.2)
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.25, top=0.97)
    save(fig, 'fig_prior')

# ------------------------------------------------------------------- Exp-5: by clique size and answer size ----
def rep_profiles():
    rows = {g: (w, x) for w, g, x in all_profiles() if g in REP}
    return [(g, rows[g][1]) for g in REP if g in rows]

def fig_profile():
    """Own-level latency by clique size (stratified workload) and by answer size (deciles), on the representative pair."""
    rows = rep_profiles()
    if not rows: return
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(TEXT_WIDTH, 1.5))
    for g, x in rows:
        bs = x['by_size']; mk = MARKS[g]
        a.plot([e['s'] for e in bs], [e['locate_ns'] for e in bs], marker=mk, ms=2.8, **OURS_KW, label=g)
        a.plot([e['s'] for e in bs], [e['st_locate_ns'] for e in bs], marker=mk, ms=2.8, **BASE_KW)
        b.plot([e['s'] for e in bs], [e['list_ns'] for e in bs], marker=mk, ms=2.8, **OURS_KW)
        b.plot([e['s'] for e in bs], [e['st_list_ns'] for e in bs], marker=mk, ms=2.8, **BASE_KW)
        dc = x['own_by_output']
        c.plot([e['output'] for e in dc], [e['list_ns'] for e in dc], marker=mk, ms=2.8, **OURS_KW)
        c.plot([e['output'] for e in dc], [e['st_list_ns'] for e in dc], marker=mk, ms=2.8, **BASE_KW)
    a.set_xscale('log'); a.set_xlabel('clique size $s$'); a.set_ylabel('locating time (ns)'); a.set_ylim(0, None)
    a.legend(handles=[Line2D([], [], marker=MARKS[g], color='0.25', ls='', markersize=3.4, label=g) for g, x in rows], loc='upper right')
    b.set_xscale('log'); b.set_yscale('log'); b.set_xlabel('clique size $s$'); b.set_ylabel('listing time (ns)')
    legend_pair(b, loc='upper right')
    c.set_xscale('log'); c.set_yscale('log'); c.set_xlabel('answer size (vertices, decile mean)'); c.set_ylabel('listing time (ns)')
    for ax, t in zip((a, b, c), ('(a)', '(b)', '(c)')): ax.set_title(t, loc='left', fontsize=8, pad=3)
    fig.subplots_adjust(left=0.07, right=0.995, bottom=0.24, top=0.89, wspace=0.42)
    save(fig, 'fig_profile')

# --------------------------------------------------------------------------------- Exp-7: scalability ----
def prior_total(where, g):
    """CND summed build and peel time (s) over all sizes for (machine, graph), from the prior records."""
    d = {'laptop': 'prior', 'tods2': 'prior/tods2', 'tods1': 'prior/tods1'}[where]
    p = EV / d / f'prior_original_{g}.json'
    if not p.exists(): p = EV / d / f"prior_original_{ {v: k for k, v in NAMES.items()}.get(g, g) }.json"
    if not p.exists(): return None
    og = json.loads(p.read_text()); return og['total_inproc_ms'] / 1000 if 'total_inproc_ms' in og else None

def fig_scale(tag='scale_tods2', full='tods2.json'):
    """Vertex-induced samples at 20 to 100 percent: bytes (index, S trees), build time (index, CND all sizes),
    own-level locating and listing time (index, S trees)."""
    rec = load(f'{tag}.json'); fullrec = load(full)
    if rec is None or fullrec is None: return
    series = {}
    for r in rec['runs']:
        if 'result' not in r: continue
        stem = Path(r['graph']).stem; g, p = stem.rsplit('_p', 1); series.setdefault(g, {})[int(p)] = r['result']
    for r in fullrec['runs']:
        if 'result' in r and Path(r['graph']).stem in series: series[Path(r['graph']).stem][100] = r['result']
    prof = {g: x for w, g, x in all_profiles(samples=True)}; prof.update({g: x for w, g, x in all_profiles() if w == 'tods2'})
    fig, axes = plt.subplots(1, 4, figsize=(TEXT_WIDTH, 1.45))
    for g, pts in sorted(series.items()):
        ps = sorted(pts); xs = [pts[p]['n'] / 1e6 for p in ps]; mk = MARKS.get(g, 'o')
        axes[0].plot(xs, [pts[p]['bytes_total'] / 1048576 for p in ps], marker=mk, ms=2.8, **OURS_KW, label=g)
        axes[0].plot(xs, [pts[p]['baseline_vertex_bytes'] / 1048576 for p in ps], marker=mk, ms=2.8, **BASE_KW)
        axes[1].plot(xs, [(pts[p]['ti_ms'] + pts[p]['build_ms'] + pts[p]['compact_ms']) / 1000 for p in ps], marker=mk, ms=2.8, **OURS_KW)
        cnd = [(pts[p]['n'] / 1e6, prior_total('tods2', g if p == 100 else f'{g}_p{p}')) for p in ps]; cnd = [(x, y) for x, y in cnd if y is not None]
        if cnd: axes[1].plot([x for x, y in cnd], [y for x, y in cnd], marker=mk, ms=2.8, **BASE_KW)
        names = [g if p == 100 else f'{g}_p{p}' for p in ps]
        have = [(pts[p]['n'] / 1e6, prof[nm]) for p, nm in zip(ps, names) if nm in prof]
        if have:
            axes[2].plot([x for x, y in have], [y['regimes']['own']['locate_ns'] for x, y in have], marker=mk, ms=2.8, **OURS_KW)
            axes[2].plot([x for x, y in have], [y['regimes']['own']['st_locate_ns'] for x, y in have], marker=mk, ms=2.8, **BASE_KW)
            axes[3].plot([x for x, y in have], [y['regimes']['own']['list_ns'] / y['regimes']['own']['output'] for x, y in have], marker=mk, ms=2.8, **OURS_KW)
            axes[3].plot([x for x, y in have], [y['regimes']['own']['st_list_ns'] / y['regimes']['own']['output'] for x, y in have], marker=mk, ms=2.8, **BASE_KW)
    for ax, lab, t in zip(axes, ('index size (MB)', 'build time, all sizes (s)', 'locating time (ns)', 'listing time per vertex (ns)'), ('(a)', '(b)', '(c)', '(d)')):
        ax.set_xlabel('vertices (millions)'); ax.set_ylabel(lab); ax.set_xlim(0, None)
        ax.set_title(t, loc='left', fontsize=8, pad=3)
    axes[0].set_yscale('log'); axes[1].set_yscale('log'); axes[2].set_ylim(0, None); axes[3].set_ylim(0, None)
    axes[0].legend(handles=[Line2D([], [], marker=MARKS.get(g, 'o'), color='0.25', ls='', markersize=3.4, label=g) for g in sorted(series)], loc='lower right')
    legend_pair(axes[3], base_label='STrees / CND', loc='lower right')
    fig.subplots_adjust(left=0.065, right=0.995, bottom=0.26, top=0.89, wspace=0.5)
    save(fig, 'fig_scale')

if __name__ == '__main__':
    fig_regimes(); fig_profile(); fig_scale(); fig_prior()
    print('figures written to', OUT)
