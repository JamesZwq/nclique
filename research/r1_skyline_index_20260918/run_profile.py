#!/usr/bin/env python3
"""Query-latency profiles (query_profile) of every index file given, one at a time, one thread.  Writes <tag>.json (a run
per index file with the program's JSON) and <tag>-logs/; refuses to overwrite.  Usage: run_profile.py --tag <tag> <index.cx> ..."""
import datetime, json, os, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]

def main():
    args = sys.argv[1:]; assert args[0] == '--tag', __doc__
    tag, files = args[1], args[2:]; out = HERE / f'{tag}.json'; logs = HERE / f'{tag}-logs'
    if out.exists() or logs.exists(): raise SystemExit('refusing to overwrite evidence')
    logs.mkdir(); build = HERE / 'build'
    for command in ([['cmake', '-S', str(HERE), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release', '-DSANITIZE=OFF'], ['cmake', '--build', str(build), '-j', '12', '--target', 'query_profile']]):
        subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
    rec = {'tag': tag, 'host': os.uname().nodename, 'started': datetime.datetime.now().astimezone().isoformat(), 'runs': []}
    for f in files:
        name = Path(f).stem; print('profiling', name, flush=True)
        child = subprocess.run([str(build / 'query_profile'), f], cwd=ROOT, text=True, capture_output=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
        (logs / f'{name}.log').write_text(child.stdout + child.stderr)
        if child.returncode: rec['runs'].append({'index': f, 'failed': child.returncode}); continue
        rec['runs'].append({'index': f, 'result': json.loads([x for x in child.stdout.splitlines() if x.startswith('{')][-1])})
        out.write_text(json.dumps(rec, indent=1) + '\n')
    rec['finished'] = datetime.datetime.now().astimezone().isoformat(); out.write_text(json.dumps(rec, indent=1) + '\n')

if __name__ == '__main__':
    main()
