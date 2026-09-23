#!/usr/bin/env python3
"""Per-graph table of the build-time ablation (RESULTS_FINAL 17.11): medians of ti + build + compact for the four
settings on every machine, the speedups, the phase shares of the final build and the peak memory."""
import json, statistics, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
SET = [('terminal', 'old'), ('terminal', 'fast'), ('tail', 'old'), ('tail', 'fast')]
tot = lambda x: (x['ti_ms'] + x['build_ms'] + x['compact_ms']) / 1000
print('| machine | graph | rounds | terminal+old s | terminal+fast s | tail+old s | tail+fast s | final vs before | every-vertex vs final | tree % | peel % | passes % | peak MB before -> final |')
print('|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|')
for where, f in (('laptop', 'ablation_laptop.json'), ('tods1', 'tods1_ablation.json'), ('tods2', 'tods2_ablation.json')):
    p = HERE / f
    if not p.exists(): continue
    runs = [r for r in json.loads(p.read_text())['runs'] if 'result' in r]
    graphs = list(dict.fromkeys(r['graph'] for r in runs))
    for g in graphs:
        rs = {s: [r for r in runs if r['graph'] == g and (r['solver'], r['treepass']) == s] for s in SET}
        if not all(rs.values()): continue
        med = {s: statistics.median(tot(r['result']) for r in rs[s]) for s in SET}
        fin = rs[('tail', 'fast')]; ft = med[('tail', 'fast')]
        sh = lambda k: 100 * statistics.median(r['result'][k] for r in fin) / 1000 / ft
        passes = 100 * statistics.median((r['result']['trees_ms'] + r['result']['chains_ms'] + r['result']['layout_ms'] + r['result']['compact_ms']) for r in fin) / 1000 / ft
        pk = lambda s: max(r['peak_rss_bytes'] or 0 for r in rs[s]) / 1048576
        print(f"| {where} | {Path(g).stem} | {len(fin)} | {med[SET[0]]:.3f} | {med[SET[1]]:.3f} | {med[SET[2]]:.3f} | {ft:.3f} | {med[SET[0]] / ft:.2f}x | {med[SET[1]] / ft:.2f}x | {sh('ti_ms'):.0f} | {sh('solve_ms'):.0f} | {passes:.0f} | {pk(SET[0]):,.0f} -> {pk(SET[3]):,.0f} |")
