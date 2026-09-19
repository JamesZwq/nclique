#!/usr/bin/env python3
"""Final chain index: Release and ASan/UBSan builds of chain_index_tool, both
selftests (brute force + disk round trip), then build/save/load/query on the
five graphs plus any extra graph paths given on the command line, one at a
time. Refuses to overwrite existing evidence.
Usage: run_final.py [--tag NAME] [--only] [graph ...]
  --tag NAME  write NAME.json and NAME-logs/ instead of final.json / final-logs/
  --only      run only the graphs given on the command line (skip the five defaults)"""
import datetime
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GRAPHS = ['data/ca-GrQc.edges', 'data/ca-HepPh.edges', 'data/com-dblp.edges',
          'graphs/web-Stanford.edges', 'graphs/amazon0302.edges']

def sha(path):
    path = Path(path)
    if hasattr(hashlib, 'file_digest'):
        return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''): h.update(block)
    return h.hexdigest()

TIME_FLAG = '-l' if sys.platform == 'darwin' else '-v'   # BSD time prints bytes and "real"; GNU time prints kbytes and "Elapsed"

def parse_time(text):
    """peak resident bytes and wall seconds from /usr/bin/time output (BSD -l or GNU -v)."""
    peak = wall = None
    for line in text.splitlines():
        t = line.split()
        if 'maximum resident set size' in line and t: peak = int(t[0])
        if 'Maximum resident set size' in line: peak = int(line.rsplit(':', 1)[1]) * 1024
        if len(t) >= 2 and t[1] == 'real': wall = float(t[0])
        if 'Elapsed (wall clock) time' in line:
            hms = line.rsplit(' ', 1)[1].split(':'); wall = sum(float(x) * 60 ** i for i, x in enumerate(reversed(hms)))
    return peak, wall

def run(command, log):
    child = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                           env={**os.environ, 'OMP_NUM_THREADS': '1',
                                'ASAN_OPTIONS': 'halt_on_error=1',
                                'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1'})
    log.write_text('$ ' + ' '.join(command) + '\n' + child.stdout + child.stderr)
    if child.returncode:
        raise SystemExit('failed: ' + ' '.join(command))
    return child.stdout.strip()

def main():
    args = sys.argv[1:]; evidence = 'final'; only = False
    if args[:1] == ['--tag']: evidence = args[1]; args = args[2:]
    if args[:1] == ['--only']: only = True; args = args[1:]
    graphs = ([] if only else GRAPHS) + args
    logs = HERE / f'{evidence}-logs'; cx = HERE / 'cx'
    if (HERE / f'{evidence}.json').exists() or logs.exists():
        raise SystemExit('refusing to overwrite evidence')
    logs.mkdir(); cx.mkdir(exist_ok=True)
    record = {'started': datetime.datetime.now().astimezone().isoformat(), 'commands': [],
              'sources': {str(p.relative_to(ROOT)): sha(p) for p in
                          (HERE / 'chain_index.hpp', HERE / 'chain_index_tool.cpp', HERE / 'count.cpp', HERE / 'common.hpp',
                           HERE / 'CMakeLists.txt', HERE / 'run_final.py', HERE / 'CHAINS.md', HERE / 'THEORY.md')},
              'selftests': {}, 'runs': []}
    for name, sanitize in [('build', 'OFF'), ('build-asan', 'ON')]:
        build = HERE / name
        commands = [['cmake', '-S', str(HERE), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release', f'-DSANITIZE={sanitize}'],
                    ['cmake', '--build', str(build), '-j', '12', '--target', 'chain_index_tool'],
                    [str(build / 'chain_index_tool'), '--selftest']]
        for i, command in enumerate(commands):
            record['commands'].append(command)
            out = run(command, logs / f'{name}-{i}.log')
            if i == 2:
                record['selftests'][name] = json.loads(out.splitlines()[-1])
    for graph in graphs:
        tag = Path(graph).stem
        command = ['/usr/bin/time', TIME_FLAG, str(HERE / 'build' / 'chain_index_tool'), '--bench', graph, str(cx / f'{tag}.cx')]
        record['commands'].append(command)
        output = run(command, logs / f'{tag}.log')
        line = next(x for x in output.splitlines() if x.startswith('{'))
        peak, wall = parse_time((logs / f'{tag}.log').read_text())
        record['runs'].append({'graph': graph, 'input_sha256': sha(ROOT / graph), 'result': json.loads(line),
                               'file_sha256': sha(cx / f'{tag}.cx'), 'peak_rss_bytes': peak, 'wall_s': wall, 'host': platform.node()})
    (HERE / f'{evidence}.json').write_text(json.dumps(record, indent=1) + '\n')

if __name__ == '__main__':
    main()
