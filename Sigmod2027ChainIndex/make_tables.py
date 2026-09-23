#!/usr/bin/env python3
"""Generate the paper's tables (tables/*.tex) from the evidence JSON files of
research/r1_skyline_index_20260918.  Every number in the paper's tables comes
from here; nothing is typed by hand.  Run from anywhere."""
import json, re
import statistics
from pathlib import Path
import build_records   # 2026-09-23: build time and memory from the build-time ablation (tail+fast = Section 7)

HERE = Path(__file__).resolve().parent
EV = HERE.parent / 'research' / 'r1_skyline_index_20260918'
OUT = HERE / 'tables'
OUT.mkdir(exist_ok=True)

def load(name):
    p = EV / name
    return json.loads(p.read_text()) if p.exists() else None

def fmt(n):
    return f'{int(round(n)):,}'

def mb(b):
    return f'{b / 1048576:.1f}' if b < 100 * 1048576 else f'{b / 1048576:,.0f}'

def tex_escape(s):
    return s.replace('\\', '/').replace('&', '\\&').replace('%', '\\%').replace('#', '\\#').replace('_', '\\_').replace('$', '\\$')

FAMILY = {'ca-': 'collab', 'dblp': 'collab', 'cit-': 'citation', 'web-': 'web', 'amazon': 'product', 'com-amazon': 'product',
          'email': 'comm.', 'loc-': 'social', 'soc-': 'social', 'com-youtube': 'social', 'wiki': 'comm.', 'tech-': 'internet'}
MACHINE = {'laptop': 'laptop', 'tods1': 'server 1', 'tods2': 'server 2'}
NAMES = {'soc-pokec-relationships': 'soc-pokec', 'com-amazon.ungraph': 'com-amazon'}

def family(g):
    if g in ('com-dblp', 'dblp-core30', 'ca-dblp-2012', 'dblp-coauthor', 'ca-coauthors-dblp'): return 'collab'
    for k, v in FAMILY.items():
        if g.startswith(k) or k in g: return v
    return 'other'

def rows_by_graph():
    """All recorded runs, keyed by distinct graph (n, m), then by machine."""
    rows = {}
    for where, f in [('tods1', 'tods1.json'), ('tods2', 'tods2.json'), ('laptop', 'final.json'), ('laptop', 'more.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' not in r: continue
            key = (r['result']['n'], r['result']['m'])
            rows.setdefault(key, {}).setdefault(where, (NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem), r))
    return rows

def merged_rows(prefer=('tods1', 'tods2', 'laptop')):
    """One row per distinct graph (by n, m), taking the first machine of `prefer` that has it."""
    best = []
    for key, by in rows_by_graph().items():
        where = next(w for w in prefer if w in by); g, r = by[where]; best.append((where, g, r))
    return sorted(best, key=lambda t: t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'])

def table_size():
    """One row per graph: chains, values stored, bytes of the index and of STrees, build time and peak memory (build table folded in)."""
    rows = merged_rows(); B = build_records.builds()
    lines = [r'\begin{tabular}{@{}llrrrrrrrrrr@{}}', r'\toprule',
             r'Graph & Type & $n$ & $s_{\max}$ & Chains & $n/$chains & Values stored & \chainidx (MB) & \strees (MB) & Ratio & Build (s) & Memory (GB) \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']; b = B.get((where, g))
        if b is None:   # no build-time ablation record yet: the older record, flagged
            print(f'PROVISIONAL build numbers for {g} on {where} (older record)'); b = {'total_s': (x['ti_ms'] + x['build_ms'] + x['compact_ms']) / 1000, 'peak_bytes': r.get('peak_rss_bytes')}
        total = b['total_s']; peak = b['peak_bytes']
        lines.append(f"{tex_escape(g)} & {family(g)} & {fmt(x['n'])} & {x['s_max']} & {fmt(x['chains'])} & {x['n']/x['chains']:.1f} & {100*x['vertex_residue_cells']/x['vertex_pairs']:.1f}\\% & {mb(x['bytes_total'])} & {mb(x['baseline_vertex_bytes'])} & {x['baseline_vertex_bytes']/x['bytes_total']:.1f}$\\times$ & {total:.2f} & {(peak/1073741824):.2f} \\\\" if peak else
                     f"{tex_escape(g)} & {family(g)} & {fmt(x['n'])} & {x['s_max']} & {fmt(x['chains'])} & {x['n']/x['chains']:.1f} & {100*x['vertex_residue_cells']/x['vertex_pairs']:.1f}\\% & {mb(x['bytes_total'])} & {mb(x['baseline_vertex_bytes'])} & {x['baseline_vertex_bytes']/x['bytes_total']:.1f}$\\times$ & {total:.2f} & -- \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'size.tex').write_text('\n'.join(lines) + '\n')
    rat = [t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'] for t in rows]
    with_perm = [t[2]['result']['baseline_vertex_bytes'] / (t[2]['result']['bytes_total'] + 4 * t[2]['result']['n']) for t in rows]
    (OUT / 'size_stats.tex').write_text(f"\\newcommand{{\\numgraphs}}{{{len(rows)}}}\n\\newcommand{{\\ratiomin}}{{{min(rat):.2f}}}\n\\newcommand{{\\ratiomedian}}{{{statistics.median(rat):.1f}}}\n\\newcommand{{\\ratiomax}}{{{max(rat):.1f}}}\n\\newcommand{{\\permmin}}{{{min(with_perm):.2f}}}\n\\newcommand{{\\permmedian}}{{{statistics.median(with_perm):.1f}}}\n\\newcommand{{\\permmax}}{{{max(with_perm):.1f}}}\n")

def strees_latency():
    """S trees latency keyed by (machine, graph): the in-process baseline of query_profile (same decomposition, same queries,
    parent-pointer climb, one memory copy per answer), from profile_<machine>.json.  Keys mirror the index's fields:
    own/half/root list and locate times in ns."""
    out = {}
    for where in ('laptop', 'tods1', 'tods2'):
        rec = load(f'profile_{where}.json')
        if rec is None: continue
        for r in rec['runs']:
            if 'result' not in r: continue
            g = NAMES.get(Path(r['index']).stem, Path(r['index']).stem)
            if '_p' in g: continue
            R = r['result']['regimes']
            out[(where, g)] = {f'{reg}_base_ns': R[reg]['st_list_ns'] for reg in R} | {f'{reg}_locate_ns': R[reg]['st_locate_ns'] for reg in R} | {f'{reg}_ours_ns': R[reg]['list_ns'] for reg in R} | {f'{reg}_ours_locate_ns': R[reg]['locate_ns'] for reg in R}
    return out

def table_queries():
    """One row per graph.  Community columns come from the profile record of the machine (index and S trees measured in one
    process on the same queries); the value query from the run_final record of the same machine."""
    st = strees_latency(); rows = []
    for key, by in rows_by_graph().items():
        with_base = [w for w in ('tods2', 'tods1', 'laptop') if w in by and (w, by[w][0]) in st]
        where = with_base[0] if with_base else next(w for w in ('tods1', 'tods2', 'laptop') if w in by); g, r = by[where]; rows.append((where, g, r))
    rows.sort(key=lambda t: t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'])
    prof = {}
    for where in ('laptop', 'tods1', 'tods2'):
        rec = load(f'profile_{where}.json')
        if rec is None: continue
        for r in rec['runs']:
            if 'result' in r: prof[(where, NAMES.get(Path(r['index']).stem, Path(r['index']).stem))] = r['result']['regimes']['own']
    lines = [r'\begin{tabular}{@{}llrrrrrrrr@{}}', r'\toprule',
             r'Graph & Machine & Locate (ns) & \strees locate (ns) & Ranges & Vertices & List (ns) & \strees list (ns) & List (ns/vertex) & Value (ns) \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']; o = prof.get((where, g))
        if o: lines.append(f"{tex_escape(g)} & {MACHINE[where]} & {o['locate_ns']:.0f} & {o['st_locate_ns']:.0f} & {o['ranges']:.0f} & {fmt(o['output'])} & {o['list_ns']:,.0f} & {o['st_list_ns']:,.0f} & {o['list_ns']/o['output']:.2f} & {x['value_ns']:.0f} \\\\")
        else: lines.append(f"{tex_escape(g)} & {MACHINE[where]} & {x['ptr_own_ns']:.0f} & -- & {x['own_ranges']:.0f} & {fmt(x['own_output'])} & {x['explicit_own_ns']:,.0f} & -- & {x['explicit_own_ns']/x['own_output']:.2f} & {x['value_ns']:.0f} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'queries.tex').write_text('\n'.join(lines) + '\n')
    # macros for the text, from the same rows as the table (profile numbers where a profile exists)
    all_rows = merged_rows(); x = lambda k: [r['result'][k] for w, g, r in all_rows]
    picked = [(where, g, prof.get((where, g)), r['result']) for where, g, r in rows]
    loc = [o['locate_ns'] if o else r['ptr_own_ns'] for w, g, o, r in picked]
    outs = [(o['list_ns'], o['output']) if o else (r['explicit_own_ns'], r['own_output']) for w, g, o, r in picked]
    big = [l / n for l, n in outs if n >= 1000]
    macros = {'locmin': f"{min(loc):.0f}", 'locmax': f"{max(loc):.0f}", 'locmedian': f"{statistics.median(loc):.0f}",
              'pervmin': f"{min(big):.2f}", 'pervmax': f"{max(big):.2f}",
              'valmin': f"{min(x('value_ns')):.0f}", 'valmax': f"{max(x('value_ns')):.0f}", 'valmedian': f"{statistics.median(x('value_ns')):.0f}",
              'storedmin': f"{min(100 * r['result']['vertex_residue_cells'] / r['result']['vertex_pairs'] for w, g, r in all_rows):.1f}",
              'storedmax': f"{max(100 * r['result']['vertex_residue_cells'] / r['result']['vertex_pairs'] for w, g, r in all_rows):.1f}",
              'storedmedian': f"{statistics.median(100 * r['result']['vertex_residue_cells'] / r['result']['vertex_pairs'] for w, g, r in all_rows):.0f}"}
    # levels, from the profiles (index and S trees in one process)
    half = [s['half_ours_locate_ns'] for s in st.values()]; root = [s['root_ours_locate_ns'] for s in st.values()]; own = [s['own_ours_locate_ns'] for s in st.values()]
    macros.update({'ownmedian': f"{statistics.median(own):.0f}", 'halfmedian': f"{statistics.median(half):.0f}", 'halfmax': f"{max(half):.0f}", 'rootmedian': f"{statistics.median(root):.0f}"})
    # index against S trees, both from the profiles
    faster = total = 0; ratios = []; own_loc = []; deep_loc = []
    for (where, g), s in st.items():
        for reg in ('own', 'half', 'root'):
            ratio = s[f'{reg}_base_ns'] / s[f'{reg}_ours_ns']; ratios.append(ratio); total += 1; faster += ratio > 1
        own_loc.append(s['own_locate_ns'] / s['own_ours_locate_ns']); deep_loc.append(max(s['half_locate_ns'] / s['half_ours_locate_ns'], s['root_locate_ns'] / s['root_ours_locate_ns']))
    macros.update({'stpoints': str(total), 'stfaster': str(faster), 'stpairs': str(len(st)), 'stgraphs': str(len({g for (w, g) in st})),
                   'stmin': f"{min(ratios):.1f}" if ratios else '?', 'stmax': f"{max(ratios):.1f}" if ratios else '?', 'stmedian': f"{statistics.median(ratios):.2f}" if ratios else '?',
                   'stmininv': f"{1/max(ratios):.2f}" if ratios else '?', 'stmaxinv': f"{1/min(ratios):.1f}" if ratios else '?', 'stmedianinv': f"{1/statistics.median(ratios):.2f}" if ratios else '?',
                   'stlocown': f"{statistics.median(own_loc):.1f}" if own_loc else '?', 'stlocdeepmax': f"{max(deep_loc):.0f}" if deep_loc else '?', 'stlocdeepmedian': f"{statistics.median(deep_loc):.1f}" if deep_loc else '?'})
    # Exp-5: the representative pair of the by-size and by-answer-size profiles (figure fig_profile)
    REP = ['web-BerkStan', 'cit-Patents']; prof = {}
    for where in ('tods2', 'tods1', 'laptop'):
        rec = load(f'profile_{where}.json')
        for r in (rec or {'runs': []})['runs']:
            g = NAMES.get(Path(r['index']).stem, Path(r['index']).stem)
            if 'result' in r and g in REP: prof.setdefault(g, r['result'])
    if len(prof) == len(REP):
        bs = [e for g in REP for e in prof[g]['by_size']]
        small = [e['st_list_ns'] / e['list_ns'] for e in bs if e['s'] <= 3]; big = [e['st_list_ns'] / e['list_ns'] for e in bs if e['s'] >= 4]
        dec = [e for g in REP for e in prof[g]['own_by_output']]
        low = [e['list_ns'] / e['st_list_ns'] for e in dec if e['output'] < 1e5]; high = [abs(1 - e['st_list_ns'] / e['list_ns']) for e in dec if e['output'] >= 1e6]
        big_e = min(bs, key=lambda e: e['st_list_ns'] / e['list_ns'] if e['s'] >= 4 else 9)
        macros.update({'proflocmin': f"{min(e['locate_ns'] for e in bs):.0f}", 'proflocmax': f"{max(e['locate_ns'] for e in bs):.0f}",
                       'profstlocmin': f"{min(e['st_locate_ns'] for e in bs):.0f}", 'profstlocmax': f"{max(e['st_locate_ns'] for e in bs):.0f}",
                       'profsmallmin': f"{min(small):.2f}", 'profsmallmax': f"{max(small):.2f}", 'profsmallinvmin': f"{1/max(small):.2f}", 'profsmallinvmax': f"{1/min(small):.2f}",
                       'profbigmax': f"{1 / min(big):.0f}", 'profbigs': str(big_e['s']), 'profbigout': f"{big_e['output']:,.0f}", 'profbigranges': f"{big_e['ranges']:,.0f}",
                       'profdeclowmin': f"{min(low):.1f}", 'profdeclowmax': f"{max(low):.1f}", 'profdechigh': f"{100 * max(high):.0f}"})
    (OUT / 'query_stats.tex').write_text(''.join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in macros.items()))

def table_build():
    rows = merged_rows()
    lines = [r'\begin{tabular}{@{}llrrrrrr@{}}', r'\toprule',
             r'Graph & Machine & Clique tree (s) & Peel, all $s$ (s) & Trees and chains (s) & Total (s) & Peak RSS (GB) & Index (MB) \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']; peak = r.get('peak_rss_bytes')
        total = (x['ti_ms'] + x['build_ms'] + x['compact_ms']) / 1000
        lines.append(f"{tex_escape(g)} & {where} & {x['ti_ms']/1000:.1f} & {x['solve_ms']/1000:.1f} & {(x['trees_ms']+x['chains_ms']+x['layout_ms']+x['compact_ms'])/1000:.1f} & {total:.1f} & {(peak/1073741824):.1f} & {mb(x['bytes_total'])} \\\\" if peak else
                     f"{tex_escape(g)} & {where} & {x['ti_ms']/1000:.1f} & {x['solve_ms']/1000:.1f} & {(x['trees_ms']+x['chains_ms']+x['layout_ms']+x['compact_ms'])/1000:.1f} & {total:.1f} & -- & {mb(x['bytes_total'])} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'build.tex').write_text('\n'.join(lines) + '\n')

def table_layouts():
    """Stage-2 ablation on the five laptop graphs: per-vertex S trees, twins, chains, aligned, and the final module."""
    files = {'vertices': 'stages/index_vertices.json', 'twins': 'stages/index.json', 'chains': 'stages/index_chains.json', 'aligned': 'stages/index_aligned.json'}
    data = {k: {Path(r['graph']).stem: r['result'] for r in load(f)['runs']} for k, f in files.items()}
    final = {Path(r['graph']).stem: r['result'] for r in load('final.json')['runs']}
    graphs = ['ca-GrQc', 'ca-HepPh', 'com-dblp', 'web-Stanford', 'amazon0302']
    lines = [r'\begin{tabular}{@{}lrrrrr@{}}', r'\toprule',
             r'Graph & $S$ trees per vertex & over twin classes & over chains & chains, aligned labels & final (runs and entries) \\', r'\midrule']
    for g in graphs:
        a = data['aligned'][g]; ab = a['bytes_shared_aligned'] + a['bytes_base_nodes'] + a['bytes_base_pairs'] + a['bytes_block_d']
        lines.append(f"{tex_escape(g)} & {fmt(data['vertices'][g]['base_with_d'])} & {fmt(data['twins'][g]['base_with_d'])} & {fmt(data['chains'][g]['base_with_d'])} & {fmt(ab)} & {fmt(final[g]['bytes_total'])} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'layouts.tex').write_text('\n'.join(lines) + '\n')
    # 2026-09-22: the design study is a paragraph, not a table; ranges over the five graphs as macros
    base = {g: data['vertices'][g]['base_with_d'] for g in graphs}
    twins = [base[g] / data['twins'][g]['base_with_d'] for g in graphs]
    chains = [base[g] / data['chains'][g]['base_with_d'] for g in graphs]
    aligned = [base[g] / (data['aligned'][g]['bytes_shared_aligned'] + data['aligned'][g]['bytes_base_nodes'] + data['aligned'][g]['bytes_base_pairs'] + data['aligned'][g]['bytes_block_d']) for g in graphs]
    fin = [base[g] / final[g]['bytes_total'] for g in graphs]
    (OUT / 'layout_stats.tex').write_text(''.join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in {
        'laytwinmin': f"{min(twins):.2f}", 'laytwinmax': f"{max(twins):.2f}", 'laychainmin': f"{min(chains):.1f}", 'laychainmax': f"{max(chains):.1f}",
        'layalignmin': f"{min(aligned):.1f}", 'layalignmax': f"{max(aligned):.1f}", 'layfinalmin': f"{min(fin):.1f}", 'layfinalmax': f"{max(fin):.1f}"}.items()))

def table_prior():
    """CND (the original single-size implementation) run once per size against the chain index build, on every machine
    that has both records; the chain index build time and peak memory are from the same machine."""
    ours = {}
    for where, f in [('laptop', 'final.json'), ('laptop', 'more.json'), ('tods1', 'tods1.json'), ('tods2', 'tods2.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' in r: ours[(where, NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem))] = r
    B = build_records.builds()
    prior = []
    for where, d in [('laptop', 'prior'), ('tods2', 'prior/tods2'), ('tods1', 'prior/tods1')]:
        for p in sorted((EV / d).glob('prior_original_*.json')):
            og = json.loads(p.read_text()); g = NAMES.get(p.stem[len('prior_original_'):], p.stem[len('prior_original_'):])
            if 'total_wall_s' in og and (where, g) in ours and (where, g) in B: prior.append((where, g, og))
    lines = [r'\begin{tabular}{@{}llrrrrrr@{}}', r'\toprule',
             r'Graph & Machine & Sizes & \multicolumn{2}{c}{\cnd, all sizes (s)} & \cnd memory & \chainidx & Ratio \\',
             r' & & & total & build and peel & one size (MB) & build (s) & build and peel \\', r'\midrule']
    ratios = []
    for where, g, og in sorted(prior, key=lambda t: (t[0] != 'laptop', t[0], ours[(t[0], t[1])]['result']['n'])):
        r = ours[(where, g)]; o = r['result']; opeak = max(e['peak_rss_bytes'] or 0 for e in og['sizes'])
        build_s = B[(where, g)]['total_s']; ratio = og['total_inproc_ms'] / 1000 / build_s; ratios.append(ratio)
        peak = B[(where, g)]['peak_bytes']; opeak_s = f'{opeak/1048576:,.0f}' if opeak else '--'
        lines.append(f"{tex_escape(g)} & {MACHINE[where]} & {og['sizes_ok']} & {og['total_wall_s']:,.1f} & {og['total_inproc_ms']/1000:,.1f} & {opeak_s} & {build_s:.2f} & {ratio:.1f}$\\times$ \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'prior.tex').write_text('\n'.join(lines) + '\n')
    # 2026-09-23: the graphs at the two ends (named in Exp-1), with their numbers of sizes
    ends = sorted((og['total_inproc_ms'] / 1000 / B[(where, g)]['total_s'], g, og['sizes_ok']) for where, g, og in prior)
    (OUT / 'prior_stats.tex').write_text(f"\\newcommand{{\\priorgraphs}}{{{len(ratios)}}}\n\\newcommand{{\\priormin}}{{{min(ratios):.1f}}}\n\\newcommand{{\\priormedian}}{{{statistics.median(ratios):.0f}}}\n\\newcommand{{\\priormax}}{{{max(ratios):.0f}}}\n"
        f"\\newcommand{{\\priormingraph}}{{\\textsf{{{ends[0][1]}}}}}\n\\newcommand{{\\priorminsizes}}{{{ends[0][2]}}}\n"
        f"\\newcommand{{\\priormaxgraph}}{{\\textsf{{{ends[-1][1]}}}}}\n\\newcommand{{\\priormaxsizes}}{{{ends[-1][2]}}}\n")

def build_stats():
    """Exp-1 and the cost paragraph of Section 7 (2026-09-23), one value per Table 1 row (its preferred machine):
    the final build (tail+fast) against the build before (terminal+old) and against every vertex peeled at every size
    with the same tree pass (terminal+fast), the time shares of the final build, and the builds the text quotes."""
    fin = build_records.builds(build_records.FINAL); bef = build_records.builds(build_records.BEFORE)
    every = build_records.builds(('terminal', 'fast'))
    keys = [(where, g) for where, g, r in merged_rows() if (where, g) in fin]
    speed = sorted((bef[k]['total_s'] / fin[k]['total_s'], k) for k in keys if k in bef)
    settle = sorted((every[k]['total_s'] / fin[k]['total_s'], k) for k in keys if k in every)
    share = lambda f: sorted(f(fin[k]) / fin[k]['total_s'] for k in keys)
    tree, peel, passes = share(lambda b: b['ti_s']), share(lambda b: b['solve_s']), share(lambda b: b['trees_s'] + b['rest_s'])
    mem = sorted(fin[k]['rss_with_ti_bytes'] / fin[k]['peak_bytes'] for k in keys)   # graph + clique tree over the peak
    def quote(g):
        k = next(k for k in keys if k[1] == g); return fin[k]
    m = {'buildrows': len(keys),
         'speedmin': f"{speed[0][0]:.2f}", 'speedmedian': f"{statistics.median(x for x, _ in speed):.2f}", 'speedmax': f"{speed[-1][0]:.1f}",
         'settlemin': f"{settle[0][0]:.2f}", 'settlemedian': f"{statistics.median(x for x, _ in settle):.2f}", 'settlemax': f"{settle[-1][0]:.1f}",
         'treesharemin': f"{100*tree[0]:.0f}", 'treesharemedian': f"{100*statistics.median(tree):.0f}", 'treesharemax': f"{100*tree[-1]:.0f}",
         'peelsharemin': f"{100*peel[0]:.0f}", 'peelsharemedian': f"{100*statistics.median(peel):.0f}", 'peelsharemax': f"{100*peel[-1]:.0f}",
         'passsharemin': f"{100*passes[0]:.0f}", 'passsharemedian': f"{100*statistics.median(passes):.0f}", 'passsharemax': f"{100*passes[-1]:.0f}",
         'memsharemin': f"{100*mem[0]:.0f}", 'memsharemedian': f"{100*statistics.median(mem):.0f}"}
    for name, g in (('grqc', 'ca-GrQc'), ('dblpcoauthor', 'dblp-coauthor'), ('comdblp', 'com-dblp'), ('webuk', 'web-uk-2005'), ('berkstan', 'web-BerkStan'), ('webit', 'web-it-2004')):
        try:
            b = quote(g); m[f'build{name}'] = f"{b['total_s']:.2f}" if b['total_s'] < 10 else f"{b['total_s']:,.0f}"; m[f'mem{name}'] = f"{b['peak_bytes']/1048576:,.0f}"
            m[f'memgb{name}'] = f"{b['peak_bytes']/1073741824:.1f}"
        except StopIteration: pass
    (OUT / 'build_stats.tex').write_text(''.join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in m.items()))
    print('build stats:', m); print('   speedups', [(round(x, 2), k) for x, k in speed]); print('   settle', [(round(x, 2), k) for x, k in settle])

def table_cases():
    """Case studies: ground-truth communities (best clique size per query) and the Amazon zoom table with named examples."""
    # ground truth: com-dblp and com-amazon (com-youtube's ground truth is not cohesive; text only)
    lines = [r'\begin{tabular}{@{}lrrrrrrrr@{}}', r'\toprule',
             r'Graph & Queries & $k$-core, best $k$ & $s=3$, best $k$ & Best fixed $s$, best $k$ & Best $s$ per query, own level & Best $(s,k)$ per query & Beats $k$-core & $\mu$s per query \\', r'\midrule']
    hist_lines = []
    for g in ('com-dblp', 'com-amazon', 'com-youtube'):
        rec = load(f'case/{g}.groundtruth.json')
        if rec is None: continue
        fx = rec['fixed']; best_fixed = max(fx.items(), key=lambda kv: kv[1]['mean_ladder_all'])   # best fixed size, at its best level
        lines.append(f"{g} & {fmt(rec['queries'])} & {fx['2']['mean_ladder_all']:.3f} & {fx['3']['mean_ladder_all']:.3f} & $s={best_fixed[0]}$: {best_fixed[1]['mean_ladder_all']:.3f} & {rec['mean_best_own']:.3f} & {rec['mean_best_all']:.3f} & {100*rec['better_than_best_core']:.0f}\\% & {1e6*rec["seconds"]/rec["queries"]:.1f} \\\\")
        h = rec['best_size_hist']; tot = sum(h.values()); share = lambda a, b: 100 * sum(v for k, v in h.items() if a <= int(k) <= b) / tot
        hist_lines.append(f"{g} & {share(2,2):.0f}\\% & {share(3,3):.0f}\\% & {share(4,4):.0f}\\% & {share(5,5):.0f}\\% & {share(6,6):.0f}\\% & {share(7,999):.0f}\\% \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'case_gt.tex').write_text('\n'.join(lines) + '\n')
    (OUT / 'case_gt_hist.tex').write_text('\n'.join([r'\begin{tabular}{@{}lrrrrrr@{}}', r'\toprule', r'Graph & $s=2$ & $s=3$ & $s=4$ & $s=5$ & $s=6$ & $s\ge7$ \\', r'\midrule'] + hist_lines + [r'\bottomrule', r'\end{tabular}']) + '\n')
    # Amazon zoom
    rec = load('case/amazon-scan.json')
    if rec is not None:
        lines = [r'\begin{tabular}{@{}rrrrrrr@{}}', r'\toprule', r'$s$ & Queries & Median size & Leaf purity & Subject purity & Dominant & Pure \\', r'\midrule']
        for s, v in rec['per_size'].items():
            if v['queries'] < 100: continue
            lines.append(f"{s} & {fmt(v['queries'])} & {fmt(v['median_size'])} & {v['mean_leaf_share']:.2f} & {v['mean_subject_share']:.2f} & {v['mean_top_subject_share']:.2f} & {100*v['pure_leaf_08']:.0f}\\% \\\\")
        lines += [r'\bottomrule', r'\end{tabular}']
        (OUT / 'case_amazon.tex').write_text('\n'.join(lines) + '\n')
    # named examples: size and leaf purity per s, members of the smallest community
    rec = load('case/amazon-queries.json')
    if rec is not None:
        lines = [r'\begin{tabular}{@{}p{0.27\linewidth}p{0.26\linewidth}p{0.43\linewidth}@{}}', r'\toprule', r'Query product & Community size (leaf purity) at $s=2,3,\dots$ & Smallest community \\', r'\midrule']
        keep = ('Kind of Blue', 'Introduction to Algorithms', 'The Fellowship of the Ring', 'The Godfather')   # 2026-09-22: four rows
        for q in rec['queries']:
            if not q['title'].startswith(keep): continue
            lv = q['levels']; sizes = ', '.join(f"{L['size']:,} ({L['leaf_share']:.2f})" for L in lv)
            last = [L for L in lv if 'members' in L]
            # 2026-09-22: whole words, no bracketed tails, at most six members; the escape runs before \dots is appended
            def short(m):
                m = re.sub(r'\s*[\(\[][^\)\]]*[\)\]]', '', m).strip()
                for sep in (':', ' - '):
                    head = m.split(sep)[0].strip()
                    if len(m) > 55 and sep in m and len(head) >= 15: m = head
                return m
            if last:   # the query product itself first, then the others in stored order
                mem = last[-1]['members']; qt = q['title'] + f" ({q['group']})"
                mem = ([m for m in mem if m == qt] + [m for m in mem if m != qt])
            members = ('; '.join(tex_escape(short(m)) for m in mem[:6]) + (' \\dots' if len(mem) > 6 else '')) if last else '--'
            lines.append(f"{tex_escape(q['title'][:60])} ({q['group']}) & {sizes} & {members} \\\\")
        lines += [r'\bottomrule', r'\end{tabular}']
        (OUT / 'case_examples.tex').write_text('\n'.join(lines) + '\n')

def table_selftest():
    rec = load('final.json'); s = rec['selftests']['build']
    (OUT / 'selftest.tex').write_text(f"\\newcommand{{\\selfgraphs}}{{{fmt(s['graphs'])}}}\n\\newcommand{{\\selfcommunities}}{{{fmt(s['community_queries'])}}}\n\\newcommand{{\\selfvalues}}{{{fmt(s['value_checks'])}}}\n\\newcommand{{\\selfmembers}}{{{fmt(s['membership_checks'])}}}\n")

if __name__ == '__main__':
    table_size(); table_queries(); table_layouts(); table_prior(); build_stats(); table_cases(); table_selftest()
    print('tables written to', OUT)
