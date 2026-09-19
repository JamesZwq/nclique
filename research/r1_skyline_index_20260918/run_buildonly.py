#!/usr/bin/env python3
"""Build-only sweep: bytes, build phases, resident memory for every graph (no query passes).
Usage: run_buildonly.py [TAG] [--only graph ...]  -> TAG.json and TAG-logs/ (default buildonly); --only skips the 13 defaults."""
import datetime, hashlib, json, os, platform, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
GRAPHS = ['data/ca-GrQc.edges', 'data/ca-HepPh.edges', 'data/com-dblp.edges', 'graphs/web-Stanford.edges', 'graphs/amazon0302.edges',
          'graphs/ca-AstroPh.edges', 'graphs/ca-CondMat.edges', 'graphs/cit-HepPh.edges', 'graphs/loc-Brightkite.edges', 'graphs/soc-Epinions1.edges',
          'graphs/soc-Slashdot0902.edges', 'graphs/com-youtube.edges', 'graphs/soc-pokec.edges']
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
def main():
    args = sys.argv[1:]; tag = 'buildonly'
    if args and not args[0].startswith('--'): tag = args[0]; args = args[1:]
    graphs = GRAPHS
    if args[:1] == ['--only']: graphs = args[1:]
    out = HERE / f'{tag}.json'; logs = HERE / f'{tag}-logs'
    if out.exists() or logs.exists(): raise SystemExit('refusing to overwrite evidence')
    logs.mkdir(); binary = HERE / 'build' / 'chain_index_tool'
    rec = {'started': datetime.datetime.now().astimezone().isoformat(), 'binary_sha256': sha(binary), 'git': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
           'sources': {str(p.relative_to(ROOT)): sha(p) for p in (HERE / 'chain_index.hpp', HERE / 'chain_index_tool.cpp', HERE / 'count.cpp', ROOT / 'research/r1_terminal_20260918/terminal.hpp')}, 'runs': []}
    for g in graphs:
        cmd = ['/usr/bin/time', TIME_FLAG, str(binary), '--build-only', g]
        child = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
        (logs / (Path(g).stem + '.log')).write_text('$ ' + ' '.join(cmd) + '\n' + child.stdout + child.stderr)
        if child.returncode:
            print(g, 'FAILED:', child.stderr.strip().splitlines()[0] if child.stderr.strip() else child.returncode, flush=True)
            rec['runs'].append({'graph': g, 'input_sha256': sha(ROOT / g), 'error': child.stderr.strip().splitlines()[0] if child.stderr.strip() else str(child.returncode), 'peak_rss_bytes': parse_time(child.stderr)[0], 'wall_s': parse_time(child.stderr)[1], 'host': platform.node()}); continue
        line = next(x for x in child.stdout.splitlines() if x.startswith('{'))
        peak, wall = parse_time(child.stderr)
        rec['runs'].append({'graph': g, 'input_sha256': sha(ROOT / g), 'result': json.loads(line), 'peak_rss_bytes': peak, 'wall_s': wall, 'host': platform.node()})
        print(g, json.loads(line)['index_bytes'], 'bytes', f'peak {peak/1048576:,.0f} MB', f'wall {wall}s', flush=True)
    out.write_text(json.dumps(rec, indent=1) + '\n')
if __name__ == '__main__': main()
