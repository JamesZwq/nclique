#!/usr/bin/env python3
"""Clean timing of the solvers of tail_check --time (index preparation + all-size solve), one process per
(graph, solver, round), never two at once, under /usr/bin/time.  Each process repeats the solve R times; the rounds
alternate the solver order (A B, then B A) to cancel drift.  Usage:
    python3 bench_tail_time.py <tag> <repeats> <solver,solver,...> <graph> [<graph> ...]
Writes <tag>.json (per graph and solver: every run in ms, the median, peak RSS) and keeps every log in <tag>-logs/."""
import json, os, platform, statistics, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
tag, repeats, solvers, graphs = sys.argv[1], int(sys.argv[2]), sys.argv[3].split(','), sys.argv[4:]
out = HERE / f'{tag}.json'; logs = HERE / f'{tag}-logs'; logs.mkdir(exist_ok=True)
timer = ['/usr/bin/time', '-l'] if platform.system() == 'Darwin' else ['/usr/bin/time', '-v']
record = json.loads(out.read_text()) if out.exists() else {'host': platform.node(), 'repeats': repeats, 'graphs': {}}
for g in graphs:
    path = g if os.path.isabs(g) else str(ROOT / g)
    name = Path(g).stem
    runs = {s: [] for s in solvers}; rss = {s: 0 for s in solvers}
    for rnd, order in enumerate([solvers, solvers[::-1]]):
        for s in order:
            cmd = timer + [str(HERE / 'build' / 'tail_check'), '--time', s, path, str(repeats)]
            p = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
            (logs / f'{name}.{s}.{rnd}.log').write_text('$ ' + ' '.join(cmd) + '\n' + p.stdout + p.stderr)
            if p.returncode != 0: print(f'{name} {s} FAILED rc={p.returncode}', flush=True); continue
            runs[s] += [float(l.split(':')[1].split()[0]) for l in p.stdout.splitlines() if ' run ' in l]
            for l in p.stderr.splitlines():
                if 'maximum resident set size' in l: rss[s] = max(rss[s], int(l.split()[0]))
                if 'Maximum resident set size' in l: rss[s] = max(rss[s], int(l.split()[-1]) * 1024)
    record['graphs'][name] = {s: {'runs_ms': runs[s], 'median_ms': statistics.median(runs[s]) if runs[s] else None,
                                  'peak_rss_bytes': rss[s]} for s in solvers}
    out.write_text(json.dumps(record, indent=1))
    base = record['graphs'][name][solvers[0]]['median_ms']
    print(f'{name:24s} ' + ' '.join(f"{s}={record['graphs'][name][s]['median_ms']:9.1f}ms"
                                    f"({base / record['graphs'][name][s]['median_ms']:5.2f}x)" for s in solvers if runs[s]), flush=True)
