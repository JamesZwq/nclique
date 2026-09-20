#!/usr/bin/env python3
"""Vertex-induced random subgraphs for the scalability experiment: every vertex is kept independently with probability p
(seed 20260921; the same seed and the same vertex order give the same sample), kept vertices are relabelled in increasing
order of their old label, and the edges among them are written in the input format (first line "n m", then one edge per
line, 0-based).  Usage: sample_vertices.py <graph.edges> <out-dir> p1 [p2 ...]   -> <out-dir>/<stem>_p<percent>.edges"""
import random, sys
from pathlib import Path

def main():
    src, outdir, ps = Path(sys.argv[1]), Path(sys.argv[2]), [float(x) for x in sys.argv[3:]]
    outdir.mkdir(parents=True, exist_ok=True)
    with open(src) as f: n, m = (int(x) for x in f.readline().split())
    for p in ps:
        rng = random.Random(20260921); keep = [rng.random() < p for _ in range(n)]
        new = [-1] * n; k = 0
        for v in range(n):
            if keep[v]: new[v] = k; k += 1
        edges = []
        with open(src) as f:
            f.readline()
            for line in f:
                a, b = line.split(); a = int(a); b = int(b)
                if keep[a] and keep[b]: edges.append((new[a], new[b]))
        out = outdir / f'{src.stem}_p{int(round(p * 100))}.edges'
        with open(out, 'w') as g:
            g.write(f'{k} {len(edges)}\n'); g.writelines(f'{a} {b}\n' for a, b in edges)
        print(out, k, len(edges), flush=True)

if __name__ == '__main__':
    main()
