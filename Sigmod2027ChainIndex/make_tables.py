#!/usr/bin/env python3
"""Generate the paper's tables (tables/*.tex) from the evidence JSON files of
research/r1_skyline_index_20260918.  Every number in the paper's tables comes
from here; nothing is typed by hand.  Run from anywhere."""
import json
import statistics
from pathlib import Path

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
    return s.replace('_', '\\_')

FAMILY = {'ca-': 'collab', 'dblp': 'collab', 'cit-': 'citation', 'web-': 'web', 'amazon': 'product', 'com-amazon': 'product',
          'email': 'comm.', 'loc-': 'social', 'soc-': 'social', 'com-youtube': 'social', 'wiki': 'comm.', 'tech-': 'internet'}
NAMES = {'soc-pokec-relationships': 'soc-pokec', 'com-amazon.ungraph': 'com-amazon'}

def family(g):
    if g in ('com-dblp', 'dblp-core30', 'ca-dblp-2012', 'dblp-coauthor', 'ca-coauthors-dblp'): return 'collab'
    for k, v in FAMILY.items():
        if g.startswith(k) or k in g: return v
    return 'other'

def merged_rows():
    """One row per distinct graph (by n, m); server rows preferred; laptop rows for the rest."""
    best = {}
    for where, f in [('tods1', 'tods1.json'), ('tods2', 'tods2.json'), ('laptop', 'final.json'), ('laptop', 'more.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' not in r: continue
            key = (r['result']['n'], r['result']['m'])
            if key not in best: best[key] = (where, NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem), r)
    return sorted(best.values(), key=lambda t: t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'])

def table_size():
    rows = merged_rows()
    lines = [r'\begin{tabular}{@{}llrrrrrrr@{}}', r'\toprule',
             r'Graph & Type & $n$ & $s_{\max}$ & Chains & $n/$chains & Index (MB) & $S$ trees (MB) & Ratio \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']
        lines.append(f"{tex_escape(g)} & {family(g)} & {fmt(x['n'])} & {x['s_max']} & {fmt(x['chains'])} & {x['n']/x['chains']:.1f} & {mb(x['bytes_total'])} & {mb(x['baseline_vertex_bytes'])} & {x['baseline_vertex_bytes']/x['bytes_total']:.1f}$\\times$ \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'size.tex').write_text('\n'.join(lines) + '\n')
    rat = [t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'] for t in rows]
    (OUT / 'size_stats.tex').write_text(f"\\newcommand{{\\numgraphs}}{{{len(rows)}}}\n\\newcommand{{\\ratiomin}}{{{min(rat):.2f}}}\n\\newcommand{{\\ratiomedian}}{{{statistics.median(rat):.1f}}}\n\\newcommand{{\\ratiomax}}{{{max(rat):.1f}}}\n")

def table_queries():
    rows = merged_rows()
    lines = [r'\begin{tabular}{@{}llrrrrrrr@{}}', r'\toprule',
             r'Graph & Machine & Locate (ns) & Ranges & Copy ranges (ns) & Vertices & List (ns) & List (ns/vertex) & Value (ns) \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']
        lines.append(f"{tex_escape(g)} & {where} & {x['ptr_own_ns']:.0f} & {x['own_ranges']:.0f} & {x['range_own_ns']:,.0f} & {fmt(x['own_output'])} & {x['explicit_own_ns']:,.0f} & {x['explicit_own_ns']/x['own_output']:.2f} & {x['value_ns']:.0f} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'queries.tex').write_text('\n'.join(lines) + '\n')

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
             r'Graph & $S$ trees per vertex & over twin classes & over chains & chains, aligned labels & final (runs, widths) \\', r'\midrule']
    for g in graphs:
        a = data['aligned'][g]; ab = a['bytes_shared_aligned'] + a['bytes_base_nodes'] + a['bytes_base_pairs'] + a['bytes_block_d']
        lines.append(f"{tex_escape(g)} & {fmt(data['vertices'][g]['base_with_d'])} & {fmt(data['twins'][g]['base_with_d'])} & {fmt(data['chains'][g]['base_with_d'])} & {fmt(ab)} & {fmt(final[g]['bytes_total'])} \\\\")
    lines += [r'\midrule', r'\multicolumn{6}{@{}l}{Community listing at the own level, ns per query (explicit vertex ids)} \\']
    for g in graphs:
        a = data['aligned'][g]
        lines.append(f"{tex_escape(g)} & {fmt(data['vertices'][g]['own_base_ns'])} & {fmt(data['twins'][g]['own_base_ns'])} & {fmt(data['chains'][g]['own_base_ns'])} & {fmt(a['aligned_own_ns'])} & {fmt(final[g]['explicit_own_ns'])} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'layouts.tex').write_text('\n'.join(lines) + '\n')

def table_prior():
    """Original and optimized single-size pipelines run once per size against the chain index build, five laptop graphs."""
    ours = {}
    for f in ('final.json', 'more.json'):
        for r in load(f)['runs']: ours[Path(r['graph']).stem] = r['result']
    bo = {Path(r['graph']).stem: r for r in load('buildonly.json')['runs']}
    lines = [r'\begin{tabular}{@{}lrrrrrrrr@{}}', r'\toprule', r'Graph & Sizes & Original, all sizes (s) & Original peak RSS (MB) & Optimized (ST\_V3), all sizes (s) & Optimized hierarchies (MB) & \chainidx build (s) & Build peak RSS (MB) & \chainidx (MB) \\', r'\midrule']
    for g in ['ca-GrQc', 'ca-HepPh', 'com-dblp', 'web-Stanford', 'amazon0302']:
        v3 = load(f'prior/prior_{g}.json'); og = load(f'prior/prior_original_{g}.json')
        if v3 is None or og is None: continue
        o = ours[g]; b = bo[g]; opeak = max(e['peak_rss_bytes'] or 0 for e in og['sizes'])
        lines.append(f"{tex_escape(g)} & {og['sizes_ok']} & {og['total_wall_s']:.1f} & {opeak/1048576:,.0f} & {v3['total_wall_s']:.1f} & {v3['total_prior_bytes']/1048576:.1f} & {(o['ti_ms']+o['build_ms']+o['compact_ms'])/1000:.2f} & {b['peak_rss_bytes']/1048576:,.0f} & {o['bytes_total']/1048576:.2f} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'prior.tex').write_text('\n'.join(lines) + '\n')

def table_selftest():
    rec = load('final.json'); s = rec['selftests']['build']
    (OUT / 'selftest.tex').write_text(f"\\newcommand{{\\selfgraphs}}{{{fmt(s['graphs'])}}}\n\\newcommand{{\\selfcommunities}}{{{fmt(s['community_queries'])}}}\n\\newcommand{{\\selfvalues}}{{{fmt(s['value_checks'])}}}\n\\newcommand{{\\selfmembers}}{{{fmt(s['membership_checks'])}}}\n")

if __name__ == '__main__':
    table_size(); table_queries(); table_build(); table_layouts(); table_prior(); table_selftest()
    print('tables written to', OUT)
