#!/usr/bin/env python3
"""Prior-tool baseline: the ORIGINAL single-size r = 1 implementation of this codebase, NCliqueVertexCoreDecomposition
(the default path of `build/bin/degeneracy_cliques <graph> 1 <s> degen` with no PIVOTER_* environment variable: SDCT_Fused
build plus the tree-mutating peel, core values only), run once per clique size s = 2 .. s_max.  Records per size the
clique-tree build and peel times, the wall time of the process and the peak RSS, and sums them into the cost of
obtaining every size with the existing tool.  No optimized variant (ST_V2/ST_V3/...) is ever run.

Usage: prior_sweep.py <graph> <s_max> [tag]   -> scripts/prior_original_<tag>.json"""
import json, os, re, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent; SRC = HERE.parent; ROOT = SRC.parent
BIN = ROOT / 'build' / 'bin' / 'degeneracy_cliques'
TIME_FLAG = '-l' if sys.platform == 'darwin' else '-v'

def parse_time(text):
    peak = wall = None
    for line in text.splitlines():
        t = line.split()
        if 'maximum resident set size' in line and t: peak = int(t[0])
        if 'Maximum resident set size' in line: peak = int(line.rsplit(':', 1)[1]) * 1024
        if len(t) >= 2 and t[1] == 'real': wall = float(t[0])
        if 'Elapsed (wall clock) time' in line:
            hms = line.rsplit(' ', 1)[1].split(':'); wall = sum(float(x) * 60 ** i for i, x in enumerate(reversed(hms)))
    return peak, wall

def took(text, label):
    m = re.search(re.escape(label) + r'.*?took: ([0-9.]+) ms', text)
    return float(m.group(1)) if m else None

def main():
    graph, smax = sys.argv[1], int(sys.argv[2]); tag = sys.argv[3] if len(sys.argv) > 3 else Path(graph).stem
    out = HERE / f'prior_original_{tag}.json'
    if out.exists(): raise SystemExit(f'refusing to overwrite {out}')
    n = int(open(ROOT / graph).readline().split()[0])
    rec = {'graph': graph, 'n': n, 's_max': smax, 'binary': str(BIN), 'mode': 'original NCliqueVertexCoreDecomposition', 'started': time.strftime('%Y-%m-%dT%H:%M:%S'), 'sizes': []}
    env = {k: v for k, v in os.environ.items() if not k.startswith('PIVOTER_')}; env['OMP_NUM_THREADS'] = '1'
    for s in range(2, smax + 1):
        cmd = ['/usr/bin/time', TIME_FLAG, str(BIN), graph, '1', str(s), 'degen']
        child = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, env=env)
        text = child.stdout + child.stderr; peak, wall = parse_time(text)
        entry = {'s': s, 'rc': child.returncode, 'build_ms': took(text, 'SDCT_Fused'), 'peel_ms': took(text, 'NucleusCoreDecomposition'), 'wall_s': wall, 'peak_rss_bytes': peak}
        rec['sizes'].append(entry)
        print(f"s={s:4d} rc={child.returncode} build {entry['build_ms']} ms peel {entry['peel_ms']} ms wall {wall} s", flush=True)
        out.write_text(json.dumps(rec, indent=1) + '\n')
    ok = [e for e in rec['sizes'] if e['rc'] == 0]
    rec['total_wall_s'] = sum(e['wall_s'] or 0 for e in rec['sizes']); rec['total_inproc_ms'] = sum((e['build_ms'] or 0) + (e['peel_ms'] or 0) for e in ok)
    rec['sizes_ok'] = len(ok); out.write_text(json.dumps(rec, indent=1) + '\n')
    print('total wall', round(rec['total_wall_s'], 1), 's; in-process', round(rec['total_inproc_ms'] / 1000, 1), 's')

if __name__ == '__main__':
    main()
