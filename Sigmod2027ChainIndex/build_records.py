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
