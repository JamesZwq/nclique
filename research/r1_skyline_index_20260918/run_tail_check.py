#!/usr/bin/env python3
"""Run tail_check --graph (correctness across the five solver configurations, one timing each; the rows are hashed
inside the timed region, so use bench_tail_time.py for timing) on a list of graphs, one at a time (never two timing runs at once), under /usr/bin/time,
appending one JSON line per graph to <tag>.json and keeping every log in <tag>-logs/.  Usage:
    python3 run_tail_check.py <tag> <graph> [<graph> ...]"""
import json, os, platform, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
tag, graphs = sys.argv[1], sys.argv[2:]
out = HERE / f'{tag}.json'; logs = HERE / f'{tag}-logs'; logs.mkdir(exist_ok=True)
timer = ['/usr/bin/time', '-l'] if platform.system() == 'Darwin' else ['/usr/bin/time', '-v']
record = json.loads(out.read_text()) if out.exists() else {'host': platform.node(), 'runs': []}
for g in graphs:
    path = g if os.path.isabs(g) else str(ROOT / g)
    name = Path(g).stem
    cmd = timer + [str(HERE / 'build' / 'tail_check'), '--graph', path]
    start = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
    (logs / f'{name}.log').write_text('$ ' + ' '.join(cmd) + '\n' + p.stdout + p.stderr)
    run = {'graph': g, 'returncode': p.returncode, 'wall_s': round(time.time() - start, 3)}
    for line in p.stdout.splitlines():
        if line.startswith('{'): run['result'] = json.loads(line)
    for line in p.stderr.splitlines():
        if 'maximum resident set size' in line: run['peak_rss_bytes'] = int(line.split()[0])
        if 'Maximum resident set size' in line: run['peak_rss_bytes'] = int(line.split()[-1]) * 1024
    record['runs'].append(run); out.write_text(json.dumps(record, indent=1))
    r = run.get('result', {})
    act = max(1, r.get('active_pairs', 1))
    print(f"{name:26s} rc={p.returncode} identical={r.get('identical')} terminal={r.get('terminal_ms', 0)/1000:8.2f}s "
          f"tail={r.get('tail_ms', 0)/1000:8.2f}s residue-peel={r.get('residue_ms', 0)/1000:8.2f}s full-peel={r.get('fullpeel_ms', 0)/1000:8.2f}s "
          f"residue={100*r.get('residue_pairs', 0)/act:5.1f}% of active pairs, sizes residue/full={r.get('tail_residue_sizes')}/{r.get('tail_full_sizes')}", flush=True)
