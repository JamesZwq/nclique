#!/usr/bin/env python3
"""
Tier-2 streaming gate (spec §27 V4.4): remove k random edges from G, build
the index ONCE on the base, insert them back one at a time through a single
prototype process (maintained index/support/cores across updates, Lemma 10
end-to-end), and verify the FINAL state equals core(G) per-vertex.

Usage:
  python3 bench_dynamic_stream.py --graph graphs/com-dblp.edges --s 5 \
      --edges 50 --tmpdir <scratch>
"""
from __future__ import annotations
import argparse, os, random, subprocess, tempfile
from pathlib import Path

from bench_dynamic_locality import read_edge_file, run_v3

PROTO = os.environ.get("DYN1S_PROTO", "./build/bin/dynamic_1s_core")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--edges", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tmpdir", default=None)
    args = ap.parse_args()

    gname = Path(args.graph).stem
    workdir = args.tmpdir or tempfile.mkdtemp(prefix="dynstream_")
    Path(workdir).mkdir(parents=True, exist_ok=True)

    n, lines = read_edge_file(args.graph)
    rng = random.Random(args.seed)
    picks = set(rng.sample(range(len(lines)), args.edges))

    base_f = os.path.join(workdir, f"stream_base_{gname}_s{args.s}.edges")
    edges_f = os.path.join(workdir, f"stream_edges_{gname}_s{args.s}.txt")
    with open(base_f, "w") as bf, open(edges_f, "w") as ef:
        bf.write(f"{n} {len(lines) - len(picks)}\n")
        for i, line in enumerate(lines):
            if i in picks:
                p = line.split()
                ef.write(f"{p[0]} {p[1]}\n")
            else:
                bf.write(line + "\n")

    print(f"[{gname} s={args.s}] streaming {args.edges} inserts", flush=True)
    print("  reference core(G) ...", flush=True)
    ref = run_v3(args.graph, args.s, f"streamref_{gname}{args.s}", workdir)
    print("  base core(G - E_k) ...", flush=True)
    base = run_v3(base_f, args.s, f"streambase_{gname}{args.s}", workdir)
    core_f = os.path.join(workdir, f"stream_core_{gname}_s{args.s}.tsv")
    with open(core_f, "w") as f:
        for x, c in base.items():
            f.write(f"{x}\t{c:.0f}\n")

    proc = subprocess.run([PROTO, base_f, str(args.s), core_f, "--edges", edges_f],
                          capture_output=True, text=True, timeout=3600)
    if proc.returncode != 0:
        print(f"  PROTOTYPE FAILED rc={proc.returncode}: {proc.stderr[-400:]}")
        return
    merged = dict(base)
    stream_done = None
    for line in proc.stdout.splitlines():
        if line.startswith("CHANGED"):
            _, x, old, new = line.split()
            merged[int(x)] = float(new)
        elif line.startswith("STREAM_DONE"):
            stream_done = line
        elif line.startswith("INDEX"):
            print(f"  {line}")
    warn = [l for l in proc.stderr.splitlines() if l.startswith("WARN")]
    mism = 0
    for x in merged.keys() | ref.keys():
        if merged.get(x, 0.0) != ref.get(x, 0.0):
            mism += 1
            if mism <= 5:
                print(f"    MISMATCH v={x} got={merged.get(x,0)} want={ref.get(x,0)}")
    total_us = 0.0
    n_edges = 0
    for line in proc.stdout.splitlines():
        if line.startswith("STATS"):
            kv = dict(t.split("=") for t in line.split()[1:])
            total_us += float(kv.get("insert_us", 0))
            n_edges += 1
    print(f"  {stream_done}")
    print(f"  WARNs: {len(warn)}")
    print(f"===== STREAMING {gname} s={args.s}: edges={n_edges} "
          f"MISMATCHES={mism} (must be 0)  total_insert={total_us/1e3:.1f}ms "
          f"mean={total_us/max(1,n_edges)/1e3:.2f}ms =====")


if __name__ == "__main__":
    main()
