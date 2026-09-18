#!/usr/bin/env python3
"""Stage-2 serial driver: Release and ASan/UBSan builds, both index selftests,
the count selftest (refactor unchanged), then the five graphs one at a time.
Refuses to overwrite existing evidence. Run from anywhere; commands run at
the repository root."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
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
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else 'twins'
    assert mode in ('twins', 'chains')
    suffix = '' if mode == 'twins' else '_' + mode
    logs = HERE / ('index-logs' if mode == 'twins' else f'index-logs{suffix}')
    if (HERE / f'index{suffix}.json').exists() or logs.exists():
        raise SystemExit('refusing to overwrite evidence')
    logs.mkdir()
    record = {'mode': mode, 'started': datetime.datetime.now().astimezone().isoformat(), 'commands': [],
              'sources': {str(p.relative_to(ROOT)): sha(p) for p in
                          (HERE / 'index.cpp', HERE / 'count.cpp', HERE / 'common.hpp', HERE / 'CMakeLists.txt',
                           HERE / 'run_index.py', HERE / 'IMPLEMENTATION2.md', HERE / 'THEORY.md')},
              'selftests': {}, 'runs': []}
    for name, sanitize in [('build', 'OFF'), ('build-asan', 'ON')]:
        build = HERE / name
        commands = [['cmake', '-S', str(HERE), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release', f'-DSANITIZE={sanitize}'],
                    ['cmake', '--build', str(build), '-j', '12'],
                    [str(build / 'index'), '--selftest'],
                    [str(build / 'count'), '--selftest']]
        for i, command in enumerate(commands):
            record['commands'].append(command)
            out = run(command, logs / f'{name}-{i}.log')
            if i >= 2:
                record['selftests'][f'{name}-{Path(command[0]).name}'] = json.loads(out.splitlines()[-1])
    for graph in GRAPHS:
        tag = Path(graph).stem
        command = ['/usr/bin/time', '-l', str(HERE / 'build' / 'index'), '--graph', graph, mode]
        record['commands'].append(command)
        output = run(command, logs / f'{tag}.log')
        line = next(x for x in output.splitlines() if x.startswith('{'))
        record['runs'].append({'graph': graph, 'input_sha256': sha(ROOT / graph), 'result': json.loads(line)})
    (HERE / f'index{suffix}.json').write_text(json.dumps(record, indent=1) + '\n')

if __name__ == '__main__':
    main()
