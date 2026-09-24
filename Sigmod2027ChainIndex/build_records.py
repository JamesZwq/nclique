"""Build times and peak memory of the chain index from the build-time ablation (2026-09-23).

run_build_ablation.py ran chain_index_tool --build-only on every machine under the four settings CHAIN_SOLVER x
CHAIN_TREEPASS (terminal+old = the build of every record before 2026-09-23; tail+fast = the build described in
Section 7).  The index is byte-identical under all four, so the older records keep their sizes and latencies and only
the build columns come from here: the median over the rounds of ti + build + compact (the formula of the tables) and
the largest peak resident set of the setting's runs."""
import json, statistics
from pathlib import Path

EV = Path(__file__).resolve().parent.parent / 'research' / 'r1_skyline_index_20260918'
FILES = {'laptop': 'ablation_laptop.json', 'tods1': 'tods1_ablation.json', 'tods2': 'tods2_ablation.json'}
NAMES = {'soc-pokec-relationships': 'soc-pokec', 'com-amazon.ungraph': 'com-amazon'}
FINAL, BEFORE = ('tail', 'fast'), ('terminal', 'old')

def _total(x): return (x['ti_ms'] + x['build_ms'] + x['compact_ms']) / 1000

def builds(setting=FINAL):
    """(machine, graph) -> dict(total_s, peak_bytes, ti_s, solve_s, trees_s, rest_s, rounds) for one setting."""
    out = {}
    for where, f in FILES.items():
        p = EV / f
        if not p.exists(): continue
        by = {}
        for r in json.loads(p.read_text())['runs']:
            if 'result' not in r or (r['solver'], r['treepass']) != setting: continue
            by.setdefault(NAMES.get(Path(r['graph']).stem, Path(r['graph']).stem), []).append(r)
        for g, rs in by.items():
            med = lambda k: statistics.median(x['result'][k] for x in rs) / 1000
            out[(where, g)] = {'total_s': statistics.median(_total(x['result']) for x in rs),
                               'peak_bytes': max(x['peak_rss_bytes'] or 0 for x in rs),
                               'ti_s': med('ti_ms'), 'solve_s': med('solve_ms'), 'trees_s': med('trees_ms'),
                               'rest_s': med('chains_ms') + med('layout_ms') + med('compact_ms'), 'rounds': len(rs),
                               'rss_with_ti_bytes': statistics.median(x['result']['rss_with_ti'] for x in rs)}
    return out

# 2026-09-24 audit: the largest clique size of every graph (omega.json, read from the stored indexes).  CND was run on the
# sizes 2 .. degeneracy + 1; the sizes above the clique number have no clique, so its totals count only 2 .. omega.
OMEGA = {k: v for k, v in json.loads((EV / 'omega.json').read_text()).items() if not k.startswith('_')}

def cnd_sizes(og, g):
    """(build and peel seconds, wall seconds, number of sizes) of a CND prior record over the sizes 2 .. omega(g); checks that
    the per-size records add up to the recorded total over all sizes."""
    ok = [e for e in og['sizes'] if e['rc'] == 0]
    full = sum(e['build_ms'] + e['peel_ms'] for e in ok)
    assert abs(full - og['total_inproc_ms']) <= 1e-6 * og['total_inproc_ms'] + 1e-3, (g, full, og['total_inproc_ms'])
    keep = [e for e in ok if e['s'] <= OMEGA[g]]
    assert len(keep) == OMEGA[g] - 1, (g, len(keep), OMEGA[g])
    return sum(e['build_ms'] + e['peel_ms'] for e in keep) / 1000, sum(e['wall_s'] for e in keep), len(keep)

def one_per_pair(entries):
    """2026-09-24: tods1 and tods2 are one machine, so some graphs have two or three CND runs on it (com-amazon, com-dblp as
    com-dblp and ca-dblp-2012, web-Google, web-Stanford).  Keep one run per (machine, graph), the graph identified by (n, m),
    preferring tods1 as Table 2 does and the name Table 2 uses.  entries: (where, g, og, n, m, ...)."""
    rank = {'tods1': 0, 'tods2': 1, 'laptop': 0}
    best = {}
    for e in sorted(entries, key=lambda e: (rank[e[0]], e[1] == 'ca-dblp-2012')):
        key = ('laptop' if e[0] == 'laptop' else 'server', e[3], e[4])
        best.setdefault(key, e)
    return list(best.values())
