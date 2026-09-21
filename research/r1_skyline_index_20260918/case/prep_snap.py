#!/usr/bin/env python3
"""SNAP graphs with ground-truth communities, prepared for the chain index: the raw ungraph.txt is relabelled to
0..n-1 in increasing order of the original ids and written as <g>.edges (header "n m"); <g>.map lists the original id
of every new label; <g>.cmty holds the top-5000 communities in the new labels, one per line.
Usage: prep_snap.py <raw-dir> <g> [<g> ...]"""
import sys
from pathlib import Path

def main():
    raw = Path(sys.argv[1]); here = Path(__file__).resolve().parent
    for g in sys.argv[2:]:
        edges = []
        with open(raw / f'{g}.ungraph.txt') as f:
            for line in f:
                if line[0] == '#': continue
                a, b = line.split(); a = int(a); b = int(b)
                if a != b: edges.append((a, b) if a < b else (b, a))
        edges = sorted(set(edges)); ids = sorted({x for e in edges for x in e}); new = {x: i for i, x in enumerate(ids)}
        with open(here / f'{g}.edges', 'w') as out:
            out.write(f'{len(ids)} {len(edges)}\n'); out.writelines(f'{new[a]} {new[b]}\n' for a, b in edges)
        (here / f'{g}.map').write_text(''.join(f'{x}\n' for x in ids))
        kept = dropped = 0
        with open(raw / f'{g}.top5000.cmty.txt') as f, open(here / f'{g}.cmty', 'w') as out:
            for line in f:
                members = [new[int(x)] for x in line.split() if int(x) in new]
                if len(members) >= 2: out.write(' '.join(map(str, members)) + '\n'); kept += 1
                else: dropped += 1
        print(g, 'vertices', len(ids), 'edges', len(edges), 'communities', kept, 'dropped', dropped, flush=True)

if __name__ == '__main__':
    main()
