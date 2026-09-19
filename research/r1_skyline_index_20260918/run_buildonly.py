#!/usr/bin/env python3
"""Build-only sweep: bytes, build phases, resident memory for every graph (no query passes). Writes buildonly.json and buildonly-logs/."""
import datetime, hashlib, json, os, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
GRAPHS = ['data/ca-GrQc.edges', 'data/ca-HepPh.edges', 'data/com-dblp.edges', 'graphs/web-Stanford.edges', 'graphs/amazon0302.edges',
          'graphs/ca-AstroPh.edges', 'graphs/ca-CondMat.edges', 'graphs/cit-HepPh.edges', 'graphs/loc-Brightkite.edges', 'graphs/soc-Epinions1.edges',
          'graphs/soc-Slashdot0902.edges', 'graphs/com-youtube.edges', 'graphs/soc-pokec.edges']
def sha(p): return hashlib.file_digest(open(p, 'rb'), 'sha256').hexdigest()
def main():
    tag = sys.argv[1] if len(sys.argv) > 1 else 'buildonly'
    out = HERE / f'{tag}.json'; logs = HERE / f'{tag}-logs'
    if out.exists() or logs.exists(): raise SystemExit('refusing to overwrite evidence')
    logs.mkdir(); binary = HERE / 'build' / 'chain_index_tool'
    rec = {'started': datetime.datetime.now().astimezone().isoformat(), 'binary_sha256': sha(binary), 'git': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
           'sources': {str(p.relative_to(ROOT)): sha(p) for p in (HERE / 'chain_index.hpp', HERE / 'chain_index_tool.cpp', HERE / 'count.cpp', ROOT / 'research/r1_terminal_20260918/terminal.hpp')}, 'runs': []}
    for g in GRAPHS:
        cmd = ['/usr/bin/time', '-l', str(binary), '--build-only', g]
        child = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
        (logs / (Path(g).stem + '.log')).write_text('$ ' + ' '.join(cmd) + '\n' + child.stdout + child.stderr)
        if child.returncode: raise SystemExit('failed: ' + g)
        line = next(x for x in child.stdout.splitlines() if x.startswith('{'))
        peak = next(int(l.split()[0]) for l in child.stderr.splitlines() if 'maximum resident' in l)
        wall = next(float(l.split()[0]) for l in child.stderr.splitlines() if l.split()[1:2] == ['real'])
        rec['runs'].append({'graph': g, 'input_sha256': sha(ROOT / g), 'result': json.loads(line), 'peak_rss_bytes': peak, 'wall_s': wall})
        print(g, json.loads(line)['index_bytes'], 'bytes', f'peak {peak/1048576:,.0f} MB', f'wall {wall}s', flush=True)
    out.write_text(json.dumps(rec, indent=1) + '\n')
if __name__ == '__main__': main()
