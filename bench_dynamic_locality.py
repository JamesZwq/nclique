#!/usr/bin/env python3
"""
Kill-or-confirm experiment for DYNAMIC (1,s)-core maintenance.

Question: after inserting ONE edge, how many vertices change their
(1,s)-core number, and by how much?  Strong locality => an incremental
algorithm has room to beat full recomputation; global churn => the
direction is dead.

Method (single-edge insertion, measured backwards):
  core(G) is computed once.  For each sampled edge e in G we compute
  core(G - e); the diff  core(G) vs core(G-e)  is exactly the effect of
  inserting e into G-e.  Insertion is monotone (cores only rise), which
  doubles as a harness sanity check: any vertex whose core DROPS after
  adding an edge indicates a bug in the harness (or algorithm).

Vertex IDs: the binary dumps cores by post-degeneracy-sort internal id,
which differs between runs.  PIVOTER_DUMP_MAPPING gives new_id ->
original_id per run; all diffs are done in original-id space.

Output: one CSV row per sampled edge:
  graph,s,u,v,n_changed,max_jump,sum_jump,n_negative
plus a summary block on stdout.

Usage:
  python3 bench_dynamic_locality.py --graph graphs/com-dblp.edges --s 5 \
      --edges 300 --workers 8 --out dynamic_locality_out
"""
from __future__ import annotations
import argparse, csv, os, random, subprocess, sys, tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

BIN = "./build/bin/degeneracy_cliques"


def read_edge_file(path: str):
    """Returns (header_n, data_lines). data_lines excludes the header."""
    with open(path) as f:
        lines = f.read().splitlines()
    header = lines[0].split()
    n = int(header[0])
    return n, lines[1:]


def run_v3(graph_path: str, s: int, tag: str, workdir: str) -> dict[int, float]:
    """Run V3 on graph_path, return {original_vertex_id: core_value}."""
    core_f = os.path.join(workdir, f"core_{tag}.tsv")
    map_f = os.path.join(workdir, f"map_{tag}.tsv")
    env = os.environ.copy()
    env.update({
        "PIVOTER_RUN_ST_V3": "1",
        "PIVOTER_DUMP_CORE": core_f,
        "PIVOTER_DUMP_MAPPING": map_f,
        "OMP_NUM_THREADS": "1",
    })
    proc = subprocess.run([BIN, graph_path, "1", str(s)],
                          env=env, capture_output=True, text=True, timeout=1800)
    if proc.returncode != 0:
        raise RuntimeError(f"binary failed on {graph_path}: {proc.stderr[-500:]}")
    # new_id -> orig_id
    new2orig: dict[int, int] = {}
    with open(map_f) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a, b = line.split()
            new2orig[int(a)] = int(b)
    cores: dict[int, float] = {}
    with open(core_f) as f:
        for line in f:
            if line.startswith("#"):
                continue
            a, b = line.split()
            cores[new2orig[int(a)]] = float(b)
    os.unlink(core_f)
    os.unlink(map_f)
    return cores


def diff_cores(base: dict[int, float], plus: dict[int, float]):
    """base = core(G-e), plus = core(G).  Missing vertex => core 0.
    Returns (n_changed, max_jump, sum_jump, n_negative, changed_list)
    where changed_list = [(v, old, new), ...]."""
    changed = []
    max_jump = 0.0
    sum_jump = 0.0
    n_negative = 0
    for v in base.keys() | plus.keys():
        b = base.get(v, 0.0)
        p = plus.get(v, 0.0)
        if p != b:
            changed.append((v, b, p))
            d = p - b
            if d < 0:
                n_negative += 1
            max_jump = max(max_jump, abs(d))
            sum_jump += abs(d)
    return len(changed), max_jump, sum_jump, n_negative, changed


# ---- worker-process globals (loaded once per worker, not per job) ----
_W = {}


def _init_worker(graph_path, s, workdir, full_cores):
    n_vertices, data_lines = read_edge_file(graph_path)
    _W.update(n=n_vertices, lines=data_lines, s=s,
              workdir=workdir, full=full_cores)


def one_edge(drop_idx: int):
    lines = _W["lines"]
    tag = f"e{drop_idx}_{os.getpid()}"
    gfile = os.path.join(_W["workdir"], f"g_{tag}.edges")
    with open(gfile, "w") as f:
        f.write(f"{_W['n']} {len(lines) - 1}\n")
        for i, line in enumerate(lines):
            if i != drop_idx:
                f.write(line + "\n")
    try:
        base_cores = run_v3(gfile, _W["s"], tag, _W["workdir"])
    finally:
        os.unlink(gfile)
    parts = lines[drop_idx].split()
    u, v = int(parts[0]), int(parts[1])
    n_changed, max_jump, sum_jump, n_neg, changed = diff_cores(base_cores, _W["full"])
    return u, v, n_changed, max_jump, sum_jump, n_neg, changed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--s", type=int, required=True)
    ap.add_argument("--edges", type=int, default=300)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="dynamic_locality_out")
    ap.add_argument("--tmpdir", default=None,
                    help="scratch dir for per-edge graph copies")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    gname = Path(args.graph).stem
    workdir = args.tmpdir or tempfile.mkdtemp(prefix="dynloc_")
    Path(workdir).mkdir(parents=True, exist_ok=True)

    n_vertices, data_lines = read_edge_file(args.graph)
    print(f"[{gname} s={args.s}] n={n_vertices}, m={len(data_lines)}", flush=True)

    print("computing core(G) once ...", flush=True)
    full_cores = run_v3(args.graph, args.s, "full", workdir)
    print(f"  vertices with core>0 dump: {len(full_cores)}", flush=True)

    rng = random.Random(args.seed)
    picks = rng.sample(range(len(data_lines)), args.edges)

    csv_path = outdir / f"locality_{gname}_s{args.s}.csv"
    rows = []
    with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=_init_worker,
            initargs=(args.graph, args.s, workdir, full_cores)) as ex:
        for i, res in enumerate(ex.map(one_edge, picks)):
            rows.append(res)
            if (i + 1) % 25 == 0:
                print(f"  {i+1}/{len(picks)} edges done", flush=True)

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["graph", "s", "u", "v", "n_changed", "max_jump",
                    "sum_jump", "n_negative"])
        for (u, v, nc, mj, sj, nn, _ch) in rows:
            w.writerow([gname, args.s, u, v, nc, mj, sj, nn])
    print(f"wrote {csv_path}")

    # Per-edge changed-vertex detail (for affected-region analysis):
    # one row per (edge, changed vertex): u,v,x,old_core,new_core
    detail_path = outdir / f"detail_{gname}_s{args.s}.csv"
    with open(detail_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["u", "v", "x", "old_core", "new_core"])
        for (u, v, _nc, _mj, _sj, _nn, ch) in rows:
            for (x, old, new) in ch:
                w.writerow([u, v, x, int(old), int(new)])
    print(f"wrote {detail_path}")

    # ---- summary ----
    changed = sorted(r[2] for r in rows)
    jumps = sorted(r[3] for r in rows)
    negs = sum(r[5] for r in rows)
    zero = sum(1 for c in changed if c == 0)
    def pct(a, q): return a[min(len(a) - 1, int(q * len(a)))]
    print(f"\n===== SUMMARY  {gname} s={args.s}  ({len(rows)} single-edge inserts) =====")
    print(f"  sanity: negative-jump vertices total = {negs}  (must be 0)")
    print(f"  #changed vertices per insert:")
    print(f"    zero-change inserts: {zero}/{len(rows)} ({100.0*zero/len(rows):.1f}%)")
    print(f"    median={pct(changed,0.5)}  p90={pct(changed,0.9)}  "
          f"p99={pct(changed,0.99)}  max={changed[-1]}")
    print(f"  max core-jump per insert:")
    print(f"    median={pct(jumps,0.5)}  p90={pct(jumps,0.9)}  max={jumps[-1]}")
    frac_gt1 = sum(1 for j in jumps if j > 1) / len(jumps)
    print(f"    inserts with jump>1: {100.0*frac_gt1:.1f}%  "
          f"(k-core's +-1 lemma would put this at 0%)")


if __name__ == "__main__":
    main()
