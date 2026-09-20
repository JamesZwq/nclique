#!/usr/bin/env python3
"""S trees latency baseline on a given graph list: the stage-2 `index` program in `vertices` mode (one tree and one DFS
array per size over vertices; a community listing is one memory copy) with the fixed query workload.  Portable (GNU or
BSD time).  Writes stages/index_vertices_<tag>.json and stages/index-logs_vertices_<tag>/; refuses to overwrite.

Usage: run_strees.py --tag <tag> <graph> [<graph> ...]"""
import datetime, hashlib, json, os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]; STAGES = HERE / 'stages'
TIME_FLAG = '-l' if sys.platform == 'darwin' else '-v'

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 24), b''): h.update(chunk)
    return h.hexdigest()

def run(command, log):
    child = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
    log.write_text('$ ' + ' '.join(command) + '\n' + child.stdout + child.stderr)
    if child.returncode: raise SystemExit('failed: ' + ' '.join(command))
    return child.stdout + child.stderr

def peak_rss(text):
    for line in text.splitlines():
        t = line.split()
        if 'maximum resident set size' in line and t: return int(t[0])
        if 'Maximum resident set size' in line: return int(line.rsplit(':', 1)[1]) * 1024
    return None

def main():
    args = sys.argv[1:]; assert args[0] == '--tag', __doc__
    tag, graphs = args[1], args[2:]
    out = STAGES / f'index_vertices_{tag}.json'; logs = STAGES / f'index-logs_vertices_{tag}'
    if out.exists() or logs.exists(): raise SystemExit('refusing to overwrite evidence')
    logs.mkdir()
    build = HERE / 'build'
    record = {'mode': 'vertices', 'tag': tag, 'host': os.uname().nodename, 'started': datetime.datetime.now().astimezone().isoformat(),
              'sources': {str(p.relative_to(ROOT)): sha(p) for p in (STAGES / 'index.cpp', HERE / 'count.cpp', HERE / 'common.hpp', HERE / 'CMakeLists.txt')},
              'commands': [], 'runs': []}
    for i, command in enumerate([['cmake', '-S', str(HERE), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release', '-DSANITIZE=OFF'],
                                 ['cmake', '--build', str(build), '-j', '12', '--target', 'index'],
                                 [str(build / 'index'), '--selftest']]):
        record['commands'].append(command); text = run(command, logs / f'build-{i}.log')
        if i == 2: record['selftest'] = json.loads([x for x in text.splitlines() if x.startswith('{')][-1])
    for graph in graphs:
        name = Path(graph).stem
        command = ['/usr/bin/time', TIME_FLAG, str(build / 'index'), '--graph', graph, 'vertices']
        record['commands'].append(command); print('running', name, flush=True)
        try:
            text = run(command, logs / f'{name}.log')
        except SystemExit as e:
            record['runs'].append({'graph': graph, 'failed': str(e)}); out.write_text(json.dumps(record, indent=1) + '\n'); continue
        line = next(x for x in text.splitlines() if x.startswith('{'))
        record['runs'].append({'graph': graph, 'input_sha256': sha(ROOT / graph if not os.path.isabs(graph) else graph), 'result': json.loads(line), 'peak_rss_bytes': peak_rss(text)})
        out.write_text(json.dumps(record, indent=1) + '\n')
    record['finished'] = datetime.datetime.now().astimezone().isoformat(); out.write_text(json.dumps(record, indent=1) + '\n')

if __name__ == '__main__':
    main()
