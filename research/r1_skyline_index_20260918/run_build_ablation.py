#!/usr/bin/env python3
"""Build-time ablation of the chain index (2026-09-23).  chain_index_tool --build-only under the four settings
CHAIN_SOLVER in {terminal, tail} x CHAIN_TREEPASS in {old, fast}; terminal+old is the build of every record before
2026-09-23, tail+fast the build from then on (all four write byte-identical indexes).  Every graph, R rounds, the
setting order reversed in every other round, one process at a time (never two timing runs at once),
OMP_NUM_THREADS=1, under /usr/bin/time.  Resumable: a (graph, setting, round) already in <tag>.json is skipped.
The build time of a run is ti_ms + build_ms + compact_ms, the formula of the paper's tables.
Usage: run_build_ablation.py <tag> <rounds> <graph> [<graph> ...]   -> <tag>.json (rewritten after every run), <tag>-logs/"""
import datetime, hashlib, json, os, platform, statistics, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
SETTINGS = [('terminal', 'old'), ('terminal', 'fast'), ('tail', 'old'), ('tail', 'fast')]
if os.environ.get('ABLATION_SETTINGS'):   # e.g. "tail:old,tail:fast" for a two-setting comparison
    SETTINGS = [tuple(x.split(':')) for x in os.environ['ABLATION_SETTINGS'].split(',')]
TIME_FLAG = '-l' if sys.platform == 'darwin' else '-v'

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()

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

def total_s(x): return (x['ti_ms'] + x['build_ms'] + x['compact_ms']) / 1000

def main():
    tag, rounds, graphs = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    out = HERE / f'{tag}.json'; logs = HERE / f'{tag}-logs'; logs.mkdir(exist_ok=True)
    binary = HERE / 'build' / 'chain_index_tool'
    rec = json.loads(out.read_text()) if out.exists() else {
        'started': datetime.datetime.now().astimezone().isoformat(), 'host': platform.node(), 'binary_sha256': sha(binary),
        'git': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
        'sources': {p.name: sha(p) for p in (HERE / 'chain_index_tool.cpp', HERE / 'tail_solver.hpp', HERE / 'treepass.hpp', HERE / 'count.cpp',
                                             ROOT / 'research/r1_terminal_20260918/terminal.hpp')}, 'runs': []}
    done = {(r['graph'], r['solver'], r['treepass'], r['round']) for r in rec['runs']}
    for g in graphs:
        path = g if os.path.isabs(g) else str(ROOT / g); name = Path(g).stem
        for rnd in range(rounds):
            for solver, treepass in (SETTINGS if rnd % 2 == 0 else SETTINGS[::-1]):
                if (g, solver, treepass, rnd) in done: continue
                cmd = ['/usr/bin/time', TIME_FLAG, str(binary), '--build-only', path]
                env = {**os.environ, 'OMP_NUM_THREADS': '1', 'CHAIN_SOLVER': solver, 'CHAIN_TREEPASS': treepass}
                child = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, env=env)
                (logs / f'{name}.{solver}.{treepass}.{rnd}.log').write_text(f'$ CHAIN_SOLVER={solver} CHAIN_TREEPASS={treepass} ' + ' '.join(cmd) + '\n' + child.stdout + child.stderr)
                peak, wall = parse_time(child.stderr)
                run = {'graph': g, 'solver': solver, 'treepass': treepass, 'round': rnd, 'returncode': child.returncode, 'peak_rss_bytes': peak, 'wall_s': wall}
                if child.returncode == 0:
                    line = next(x for x in child.stdout.splitlines() if x.startswith('{')); run['result'] = json.loads(line)
                    run['input_sha256'] = sha(path)
                else: run['error'] = (child.stderr.strip().splitlines() or [str(child.returncode)])[0]
                rec['runs'].append(run); out.write_text(json.dumps(rec, indent=1) + '\n')
                print(f"{name:24s} r{rnd} {solver:8s} {treepass:4s} " + (f"total {total_s(run['result']):9.2f}s ti {run['result']['ti_ms']/1000:8.2f} solve {run['result']['solve_ms']/1000:8.2f} trees {run['result']['trees_ms']/1000:8.2f} peak {peak/1048576:8.0f} MB"
                      if child.returncode == 0 else f"FAILED {run['error']}"), flush=True)
        ok = [r for r in rec['runs'] if r['graph'] == g and 'result' in r]
        med = {f'{s}+{t}': statistics.median(total_s(r['result']) for r in ok if (r['solver'], r['treepass']) == (s, t)) for s, t in SETTINGS if any((r['solver'], r['treepass']) == (s, t) for r in ok)}
        base = med.get('terminal+old')
        print(f"== {name}: " + ' '.join(f"{k} {v:.2f}s" + (f" ({base / v:.2f}x)" if base else '') for k, v in med.items()), flush=True)

if __name__ == '__main__': main()
