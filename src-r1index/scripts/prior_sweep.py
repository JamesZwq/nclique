#!/usr/bin/env python3
"""Prior-tool baseline: the production single-size pipeline (src/, `degeneracy_cliques <graph> 1 <s> degen` with
PIVOTER_RUN_ST_V3=1 and PIVOTER_DUMP_HIER) run once per clique size s = 2 .. s_max.  Records per size the CPI build,
peel and hierarchy-build times, the wall time of the process, the peak RSS, and the number of hierarchy branches
(rows of the dumped CSV); sums them into the cost of obtaining every size's hierarchy with the existing tool.
Bytes of the prior representation: HierarchyIndexNode is 32 bytes (id, k_birth, k_death, parent, size_birth,
size_death) and owner[] is 4 bytes per vertex, per size.

Usage: prior_sweep.py <graph> <s_max> [tag]   -> scripts/prior_<tag>.json, hierarchy CSVs discarded."""
import json, os, re, subprocess, sys, tempfile, time
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
    out = HERE / f'prior_{tag}.json'
    if out.exists(): raise SystemExit(f'refusing to overwrite {out}')
    n = int(open(ROOT / graph).readline().split()[0])
    rec = {'graph': graph, 'n': n, 's_max': smax, 'binary': str(BIN), 'started': time.strftime('%Y-%m-%dT%H:%M:%S'), 'sizes': []}
    with tempfile.TemporaryDirectory() as tmp:
        for s in range(2, smax + 1):
            hier = Path(tmp) / 'hier.csv'
            env = {**os.environ, 'OMP_NUM_THREADS': '1', 'PIVOTER_RUN_ST_V3': '1', 'PIVOTER_DUMP_HIER': str(hier)}
            cmd = ['/usr/bin/time', TIME_FLAG, str(BIN), graph, '1', str(s), 'degen']
            child = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, env=env)
            text = child.stdout + child.stderr; peak, wall = parse_time(text)
            rows = sum(1 for _ in open(hier)) - 1 if hier.exists() else None
            if hier.exists(): hier.unlink()
            entry = {'s': s, 'rc': child.returncode, 'build_ms': took(text, 'ST_V3 Build'), 'peel_ms': took(text, 'ST_V3 r=1 (peel)'),
                     'hier_ms': took(text, 'hier r=1'), 'wall_s': wall, 'peak_rss_bytes': peak, 'branches': rows,
                     'prior_bytes': (rows * 32 + 4 * n) if rows is not None else None}
            rec['sizes'].append(entry)
            print(f"s={s:4d} rc={child.returncode} build {entry['build_ms']} ms peel {entry['peel_ms']} ms hier {entry['hier_ms']} ms wall {wall} s branches {rows}", flush=True)
            out.write_text(json.dumps(rec, indent=1) + '\n')
    ok = [e for e in rec['sizes'] if e['rc'] == 0 and e['branches'] is not None]
    rec['total_wall_s'] = sum(e['wall_s'] or 0 for e in rec['sizes']); rec['total_prior_bytes'] = sum(e['prior_bytes'] for e in ok)
    rec['total_inproc_ms'] = sum((e['build_ms'] or 0) + (e['peel_ms'] or 0) + (e['hier_ms'] or 0) for e in ok)
    rec['sizes_ok'] = len(ok); out.write_text(json.dumps(rec, indent=1) + '\n')
    print('total wall', round(rec['total_wall_s'], 1), 's; in-process', round(rec['total_inproc_ms'] / 1000, 1), 's; prior bytes', rec['total_prior_bytes'])

if __name__ == '__main__':
    main()
