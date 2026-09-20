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
    rows = merged_rows()
    lines = [r'\begin{tabular}{@{}llrrrrrrrr@{}}', r'\toprule',
             r'Graph & Type & $n$ & $s_{\max}$ & Chains & $n/$chains & Values stored & Index (MB) & $S$ trees (MB) & Ratio \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']
        lines.append(f"{tex_escape(g)} & {family(g)} & {fmt(x['n'])} & {x['s_max']} & {fmt(x['chains'])} & {x['n']/x['chains']:.1f} & {100*x['vertex_residue_cells']/x['vertex_pairs']:.1f}\\% & {mb(x['bytes_total'])} & {mb(x['baseline_vertex_bytes'])} & {x['baseline_vertex_bytes']/x['bytes_total']:.1f}$\\times$ \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'size.tex').write_text('\n'.join(lines) + '\n')
    rat = [t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'] for t in rows]
    (OUT / 'size_stats.tex').write_text(f"\\newcommand{{\\numgraphs}}{{{len(rows)}}}\n\\newcommand{{\\ratiomin}}{{{min(rat):.2f}}}\n\\newcommand{{\\ratiomedian}}{{{statistics.median(rat):.1f}}}\n\\newcommand{{\\ratiomax}}{{{max(rat):.1f}}}\n")

def strees_latency():
    """Own-level listing latency of S trees (stage-2 `vertices` mode) keyed by (machine, graph)."""
    out = {}
    for where, f in [('laptop', 'stages/index_vertices.json'), ('tods2', 'stages/index_vertices_tods2.json'), ('tods1', 'stages/index_vertices_tods1.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' in r: out[(where, NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem))] = r['result']
    return out

def table_queries():
    """One row per graph; the machine is the one on which the S trees latency was also measured when there is one."""
    st = strees_latency(); rows = []
    for key, by in rows_by_graph().items():
        with_base = [w for w in ('tods2', 'laptop', 'tods1') if w in by and (w, by[w][0]) in st]
        where = with_base[0] if with_base else next(w for w in ('tods1', 'tods2', 'laptop') if w in by); g, r = by[where]; rows.append((where, g, r))
    rows.sort(key=lambda t: t[2]['result']['baseline_vertex_bytes'] / t[2]['result']['bytes_total'])
    lines = [r'\begin{tabular}{@{}llrrrrrrrr@{}}', r'\toprule',
             r'Graph & Machine & Locate (ns) & Ranges & Copy ranges (ns) & Vertices & List (ns) & List (ns/vertex) & \strees list (ns) & Value (ns) \\', r'\midrule']
    for where, g, r in rows:
        x = r['result']; s = st.get((where, g))
        base = f"{s['own_base_ns']:,.0f}" if s else '--'
        lines.append(f"{tex_escape(g)} & {where} & {x['ptr_own_ns']:.0f} & {x['own_ranges']:.0f} & {x['range_own_ns']:,.0f} & {fmt(x['own_output'])} & {x['explicit_own_ns']:,.0f} & {x['explicit_own_ns']/x['own_output']:.2f} & {base} & {x['value_ns']:.0f} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'queries.tex').write_text('\n'.join(lines) + '\n')
    # macros for the text: ranges and medians over the merged rows
    all_rows = merged_rows(); x = lambda k: [r['result'][k] for w, g, r in all_rows]
    big = [r['result'] for w, g, r in all_rows if r['result']['own_output'] >= 1000]
    macros = {'locmin': f"{min(x('ptr_own_ns')):.0f}", 'locmax': f"{max(x('ptr_own_ns')):.0f}", 'locmedian': f"{statistics.median(x('ptr_own_ns')):.0f}",
              'halfmedian': f"{statistics.median(x('ptr_half_ns')):.0f}", 'halfmax': f"{max(x('ptr_half_ns')):.0f}", 'rootmedian': f"{statistics.median(x('ptr_root_ns')):.0f}",
              'pervmin': f"{min(r['explicit_own_ns'] / r['own_output'] for r in big):.2f}", 'pervmax': f"{max(r['explicit_own_ns'] / r['own_output'] for r in big):.2f}",
              'valmin': f"{min(x('value_ns')):.0f}", 'valmax': f"{max(x('value_ns')):.0f}", 'valmedian': f"{statistics.median(x('value_ns')):.0f}",
              'storedmin': f"{min(100 * r['result']['vertex_residue_cells'] / r['result']['vertex_pairs'] for w, g, r in all_rows):.1f}",
              'storedmax': f"{max(100 * r['result']['vertex_residue_cells'] / r['result']['vertex_pairs'] for w, g, r in all_rows):.1f}",
              'storedmedian': f"{statistics.median(100 * r['result']['vertex_residue_cells'] / r['result']['vertex_pairs'] for w, g, r in all_rows):.0f}"}
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
    """CND (the original single-size implementation) run once per size against the chain index build, on every machine
    that has both records; the chain index build time and peak memory are from the same machine."""
    ours = {}
    for where, f in [('laptop', 'final.json'), ('laptop', 'more.json'), ('tods1', 'tods1.json'), ('tods2', 'tods2.json')]:
        rec = load(f)
        if rec is None: continue
        for r in rec['runs']:
            if 'result' in r: ours[(where, NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem))] = r
    bo = {Path(r['graph']).stem: r for r in load('buildonly.json')['runs']}   # laptop build-only peaks
    prior = []
    for where, d in [('laptop', 'prior'), ('tods2', 'prior/tods2'), ('tods1', 'prior/tods1')]:
        for p in sorted((EV / d).glob('prior_original_*.json')):
            og = json.loads(p.read_text()); g = NAMES.get(p.stem[len('prior_original_'):], p.stem[len('prior_original_'):])
            if 'total_wall_s' in og and (where, g) in ours: prior.append((where, g, og))
    lines = [r'\begin{tabular}{@{}llrrrrrrrr@{}}', r'\toprule', r'Graph & Machine & Sizes & \cnd, all sizes (s) & of which build and peel (s) & \cnd peak RSS, one size (MB) & \chainidx build, all sizes (s) & Ratio (build and peel) & Build peak RSS (MB) & \chainidx (MB) \\', r'\midrule']
    ratios = []
    for where, g, og in sorted(prior, key=lambda t: (t[0] != 'laptop', t[0], ours[(t[0], t[1])]['result']['n'])):
        r = ours[(where, g)]; o = r['result']; opeak = max(e['peak_rss_bytes'] or 0 for e in og['sizes'])
        build_s = (o['ti_ms'] + o['build_ms'] + o['compact_ms']) / 1000; ratio = og['total_inproc_ms'] / 1000 / build_s; ratios.append(ratio)
        peak = bo[g]['peak_rss_bytes'] if where == 'laptop' and g in bo else r.get('peak_rss_bytes')
        lines.append(f"{tex_escape(g)} & {where} & {og['sizes_ok']} & {og['total_wall_s']:,.1f} & {og['total_inproc_ms']/1000:,.1f} & {opeak/1048576:,.0f} & {build_s:.2f} & {ratio:.1f}$\\times$ & {peak/1048576:,.0f} & {o['bytes_total']/1048576:.2f} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'prior.tex').write_text('\n'.join(lines) + '\n')
    (OUT / 'prior_stats.tex').write_text(f"\\newcommand{{\\priorgraphs}}{{{len(ratios)}}}\n\\newcommand{{\\priormin}}{{{min(ratios):.1f}}}\n\\newcommand{{\\priormedian}}{{{statistics.median(ratios):.0f}}}\n\\newcommand{{\\priormax}}{{{max(ratios):.0f}}}\n")

def table_compress():
    """General-purpose compression as a yardstick: xz -9e on the S trees (degeneracy labels and aligned labels) and on the index file."""
    rec = load('compress.json')
    if rec is None or not rec.get('finished'): return
    lines = [r'\begin{tabular}{@{}lrrrrrrr@{}}', r'\toprule',
             r'Graph & \strees (MB) & xz, input labels (MB) & xz, aligned labels (MB) & \chainidx (MB) & \chainidx file, xz (MB) & xz(\strees)/\chainidx & \chainidx/xz(\chainidx) \\', r'\midrule']
    for g, e in sorted(rec['graphs'].items(), key=lambda kv: kv[1]['dump']['n']):
        d = e['dump']; sd = e['strees_compressed_degeneracy']['xz9e']; sa = e['strees_compressed_aligned']['xz9e']; ix = e['index_bytes_total']; cx = e['index_compressed']['xz9e']
        lines.append(f"{tex_escape(NAMES.get(g, g))} & {d['strees_bytes']/1048576:.2f} & {sd/1048576:.2f} & {sa/1048576:.2f} & {ix/1048576:.2f} & {cx/1048576:.2f} & {sd/ix:.2f} & {e['index_file_bytes']/cx:.1f} \\\\")
    lines += [r'\bottomrule', r'\end{tabular}']
    (OUT / 'compress.tex').write_text('\n'.join(lines) + '\n')

def table_selftest():
    rec = load('final.json'); s = rec['selftests']['build']
    (OUT / 'selftest.tex').write_text(f"\\newcommand{{\\selfgraphs}}{{{fmt(s['graphs'])}}}\n\\newcommand{{\\selfcommunities}}{{{fmt(s['community_queries'])}}}\n\\newcommand{{\\selfvalues}}{{{fmt(s['value_checks'])}}}\n\\newcommand{{\\selfmembers}}{{{fmt(s['membership_checks'])}}}\n")

if __name__ == '__main__':
    table_size(); table_queries(); table_build(); table_layouts(); table_prior(); table_compress(); table_selftest()
    print('tables written to', OUT)
