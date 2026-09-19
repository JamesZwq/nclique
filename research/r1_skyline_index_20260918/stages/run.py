#!/usr/bin/env python3
"""Reproducible serial driver. It deliberately stops before graph runs on a gate failure."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent          # stages/: evidence of this stage lives here
SRC = HERE.parent                                # the CMake project (count.cpp, CMakeLists.txt, THEORY.md)
ROOT = SRC.parents[1]
GRAPHS = ['data/ca-GrQc.edges', 'data/ca-HepPh.edges', 'data/com-dblp.edges',
          'graphs/web-Stanford.edges', 'graphs/amazon0302.edges']

def sha(path):
    return hashlib.file_digest(path.open('rb'), 'sha256').hexdigest()

def run(command, log):
    child = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                           env={**os.environ, 'OMP_NUM_THREADS': '1',
                                'ASAN_OPTIONS': 'halt_on_error=1',  # detect_leaks is unsupported on macOS
                                'UBSAN_OPTIONS': 'halt_on_error=1:print_stacktrace=1'})
    log.write_text('$ ' + ' '.join(command) + '\n' + child.stdout + child.stderr)
    if child.returncode:
        raise SystemExit('failed: ' + ' '.join(command))
    return child.stdout.strip()

def main():
    logs = HERE / 'counts-logs'
    if (HERE / 'counts.json').exists() or logs.exists():
        raise SystemExit('refusing to overwrite evidence')
    logs.mkdir()
    record = {'started': datetime.datetime.now().astimezone().isoformat(), 'commands': [],
              'sources': {str(p.relative_to(ROOT)): sha(p) for p in
                          (SRC / 'count.cpp', SRC / 'CMakeLists.txt', HERE / 'run.py',
                           HERE / 'IMPLEMENTATION.md', SRC / 'THEORY.md')}, 'runs': []}
    for name, sanitize in [('build', 'OFF'), ('build-asan', 'ON')]:
        build = SRC / name
        commands = [['cmake', '-S', str(SRC), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release', f'-DSANITIZE={sanitize}'],
                    ['cmake', '--build', str(build), '-j', '12'], [str(build / 'count'), '--selftest']]
        for i, command in enumerate(commands):
            record['commands'].append(command)
            run(command, logs / f'{name}-{i}.log')
    for graph in GRAPHS:
        tag = Path(graph).stem
        command = ['/usr/bin/time', '-l', str(SRC / 'build' / 'count'), '--graph', graph]
        record['commands'].append(command)
        output = run(command, logs / f'{tag}.log')
        record['runs'].append({'graph': graph, 'input_sha256': sha(ROOT / graph), 'result': json.loads(output)})
    (HERE / 'counts.json').write_text('\n'.join(json.dumps(x) for x in record['runs']) + '\n')

if __name__ == '__main__':
    main()
