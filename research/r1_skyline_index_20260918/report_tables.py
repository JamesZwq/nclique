#!/usr/bin/env python3
"""Print the RESULTS_FINAL.md tables from the stage-2 layout evidence (index*.json) and the final module evidence (final.json)."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
LAYOUTS = [('vertices', 'stages/index_vertices.json'), ('twins', 'stages/index.json'), ('chains', 'stages/index_chains.json'), ('aligned', 'stages/index_aligned.json')]

def load(name):
    p = HERE / name
    return json.loads(p.read_text()) if p.exists() else None

def by_graph(record):
    return {Path(r['graph']).stem: r['result'] for r in record['runs'] if 'result' in r}

def fmt(x, digits=0):
    if isinstance(x, float) and digits:
        return f'{x:,.{digits}f}'
    return f'{round(x):,}'

def layouts():
    data = {k: by_graph(load(f)) for k, f in LAYOUTS}
    graphs = list(data['aligned'].keys())
    print('### Bytes with values (Block D), stage-2 layouts')
    print('| Graph | per-vertex S trees | over twins | over chains | aligned chains (index.cpp) | aligned vs per-vertex |')
    print('|---|---:|---:|---:|---:|---:|')
    for g in graphs:
        v = data['vertices'][g]['base_with_d']; t = data['twins'][g]['base_with_d']; c = data['chains'][g]['base_with_d']
        a = data['aligned'][g]; ab = a['bytes_shared_aligned'] + a['bytes_base_nodes'] + a['bytes_base_pairs'] + a['bytes_block_d']
        print(f'| {g} | {fmt(v)} | {fmt(t)} | {fmt(c)} | {fmt(ab)} | {v/ab:.2f}x |')
    print()
    print('### Community listing, own level (k = kappa_s(v)), ns per query, explicit ids unless noted')
    print('| Graph | output vertices | per-vertex S trees (memcpy) | over twins | over chains | aligned explicit | aligned ranges only |')
    print('|---|---:|---:|---:|---:|---:|---:|')
    for g in graphs:
        a = data['aligned'][g]
        print(f"| {g} | {fmt(a['own_output'])} | {fmt(data['vertices'][g]['own_base_ns'])} | {fmt(data['twins'][g]['own_base_ns'])} | {fmt(data['chains'][g]['own_base_ns'])} | {fmt(a['aligned_own_ns'])} | {fmt(a['aligned_range_own_ns'])} |")
    print()
    print('### Value queries, ns per query (stage-2 layouts)')
    print('| Graph | per-vertex | chains | aligned |')
    print('|---|---:|---:|---:|')
    for g in graphs:
        a = data['aligned'][g]
        print(f"| {g} | {data['vertices'][g]['value_ns']:.1f} | {data['chains'][g]['value_ns']:.1f} | {a['value_ns']:.1f} |")
    print()

def final():
    rec = load('final.json')
    if rec is None:
        print('(final.json not present yet)'); return
    data = by_graph(rec); vert = by_graph(load('stages/index_vertices.json'))
    more = load('more.json')
    if more is not None: data.update(by_graph(more))   # the four graphs added on 2026-09-19 (ca-HepTh, email-Eu-core, com-amazon, dblp-coauthor)
    def vbytes(g, r): return r.get('baseline_vertex_bytes', vert[g]['base_with_d'] if g in vert else 0)
    print('### Final module: size')
    print('| Graph | n | s_max | W bits | chains | n / chains | canonical nodes | chains / nodes | (chain,s) pairs | runs | map B | chains B | layers B | total B | file B | perm B | build form total B | per-vertex S trees B | ratio | ratio with perm |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for g, r in data.items():
        v = vbytes(g, r)
        print(f"| {g} | {fmt(r['n'])} | {r['s_max']} | {r['count_bits']} | {fmt(r['chains'])} | {r['n']/r['chains']:.1f} | {fmt(r['canonical_nodes'])} | {r['chains']/r['canonical_nodes']:.2f} | {fmt(r['pairs_total'])} | {fmt(r['runs_total'])} | {fmt(r['bytes_map'])} | {fmt(r['bytes_chains'])} | {fmt(r['bytes_layers'])} | {fmt(r['bytes_total'])} | {fmt(r['file_bytes'])} | {fmt(r['perm_bytes'])} | {fmt(r['slice_bytes_total'])} | {fmt(v)} | {v/r['bytes_total']:.2f}x | {v/(r['bytes_total']+r['perm_bytes']):.2f}x |")
    print()
    def wall_rss(g):
        p = HERE / 'final-logs' / f'{g}.log'
        if not p.exists(): p = HERE / 'more-logs' / f'{g}.log'
        if not p.exists(): return '-', '-'
        wall = rss = '-'
        for line in p.read_text().splitlines():
            t = line.split()
            if len(t) >= 2 and t[1] == 'real': wall = t[0]
            if 'Elapsed (wall clock) time' in line:
                hms = line.rsplit(' ', 1)[1].split(':'); wall = f"{sum(float(x) * 60 ** i for i, x in enumerate(reversed(hms))):.2f}"
            if 'maximum resident set size' in line: rss = f'{int(t[0])/1048576:,.0f}'
            if 'Maximum resident set size' in line: rss = f'{int(line.rsplit(":", 1)[1])/1024:,.0f}'
        return wall, rss
    print('### Final module: build, save, load (ms, single thread); whole-process wall time (s) and peak RSS (MB)')
    print('| Graph | solve (all-size peel) | trees | chains + labels | layout | build total | compact | save | load | process wall s | peak RSS MB |')
    print('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for g, r in data.items():
        wall, rss = wall_rss(g)
        print(f"| {g} | {fmt(r['solve_ms'])} | {fmt(r['trees_ms'])} | {fmt(r['chains_ms'])} | {fmt(r['layout_ms'])} | {fmt(r['build_ms'])} | {r['compact_ms']:.2f} | {r['save_ms']:.1f} | {r['load_ms']:.1f} | {wall} | {rss} |")
    print()
    print('### Final module: community queries, ns per query (compact form, loaded from disk)')
    print('| Graph | regime | output vertices | ranges | locate (pointer) | ranges copied | explicit ids | per-vertex S trees memcpy | build form ranges | build form explicit |')
    print('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for g, r in data.items():
        for reg, key in [('own', 'own'), ('half', 'half'), ('root', 'root')]:
            memcpy = fmt(vert[g][f'{key}_base_ns']) if g in vert else '-'
            print(f"| {g} | {reg} | {fmt(r[f'{key}_output'])} | {r[f'{key}_ranges']:.1f} | {r[f'ptr_{key}_ns']:.1f} | {fmt(r[f'range_{key}_ns'])} | {fmt(r[f'explicit_{key}_ns'])} | {memcpy} | {fmt(r[f'slice_range_{key}_ns'])} | {fmt(r[f'slice_explicit_{key}_ns'])} |")
    print()
    print('### Final module: value queries (ns per query)')
    print('| Graph | value | max tree depth |')
    print('|---|---:|---:|')
    for g, r in data.items():
        print(f"| {g} | {r['value_ns']:.1f} | {r['max_depth']} |")
    print()
    if all('full_bytes_total' in r for r in data.values()):
        print('### Top encoding ablation, same process: compact form with tops as T against packed tops (ns per query)')
        print('| Graph | bytes T tops | bytes packed | climb own T / packed | climb half T / packed | climb root T / packed | locate own T / packed | value T / packed |')
        print('|---|---:|---:|---:|---:|---:|---:|---:|')
        for g, r in data.items():
            print(f"| {g} | {fmt(r['full_bytes_total'])} | {fmt(r['bytes_total'])} | {r['full_climb_own_ns']:.1f} / {r['climb_own_ns']:.1f} | {r['full_climb_half_ns']:.1f} / {r['climb_half_ns']:.1f} | {r['full_climb_root_ns']:.1f} / {r['climb_root_ns']:.1f} | {r['full_ptr_own_ns']:.1f} / {r['ptr_own_ns']:.1f} | {r['full_value_ns']:.1f} / {r['value_ns']:.1f} |")
        print()
        print('### Climb only (own node lookup + climb), ns per query: build form (T tops, chain-id arrays) against the loaded packed index')
        print('| Graph | own build / packed | half build / packed | root build / packed |')
        print('|---|---:|---:|---:|')
        for g, r in data.items():
            print(f"| {g} | {r['slice_climb_own_ns']:.1f} / {r['climb_own_ns']:.1f} | {r['slice_climb_half_ns']:.1f} / {r['climb_half_ns']:.1f} | {r['slice_climb_root_ns']:.1f} / {r['climb_root_ns']:.1f} |")
        print()
    print('### Explicit listing cost per output vertex (ns)')
    print('| Graph | final module (own) | per-vertex S trees (own) | final module (root) | per-vertex S trees (root) |')
    print('|---|---:|---:|---:|---:|')
    for g, r in data.items():
        if g not in vert: print(f"| {g} | {r['explicit_own_ns']/r['own_output']:.3f} | - | {r['explicit_root_ns']/r['root_output']:.3f} | - |"); continue
        print(f"| {g} | {r['explicit_own_ns']/r['own_output']:.3f} | {vert[g]['own_base_ns']/vert[g]['own_output']:.3f} | {r['explicit_root_ns']/r['root_output']:.3f} | {vert[g]['root_base_ns']/vert[g]['root_output']:.3f} |")
    print()
    print('selftests (final.json):', json.dumps(rec['selftests']))
    if more is not None: print('selftests (more.json):', json.dumps(more['selftests']))
    ratios = sorted(vbytes(g, r) / r['bytes_total'] for g, r in data.items())
    print(f"byte ratio over {len(ratios)} graphs: min {ratios[0]:.2f}x, median {ratios[len(ratios)//2]:.2f}x, max {ratios[-1]:.2f}x")

if __name__ == '__main__':
    layouts(); final()
