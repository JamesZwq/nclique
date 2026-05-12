#!/usr/bin/env python3
"""
Scalability benchmark: vertex-induced subsampling on a fixed graph.

For each ratio in `ratios`, we keep `round(ratio * n)` vertices selected
uniformly at random (nested: 20% ⊂ 40% ⊂ 60% ⊂ 80% ⊂ 100% by construction,
because the same seeded RNG yields the same shuffled vertex order every
call), then write the induced subgraph (only edges with both endpoints
kept). For each subsampled graph we run Ours_ST + REF_R1 across a sweep
of s values and report wall-clock time and peak RSS.

Why vertex-induced rather than random-edge: SPIN★ / CND runtime scales
with Σ (encoded s-clique count), and random edge removal drops Σ
super-linearly because each missing edge can destroy many cliques. A
vertex-induced subsample preserves clique structure inside the kept
vertex set; the # of s-cliques scales smoothly with the kept vertex
count, which gives clean monotone time / memory curves.

Usage:
    python3 bench_r1_scalability.py \
        --bin ./build/bin/degeneracy_cliques \
        --graph graphs/web-Google.edges \
        --ratios 0.2 0.4 0.6 0.8 1.0 \
        --s-list 2 3 5 7 9 11 13 15 \
        --runs 1 \
        --out paper_data/scalability.csv
"""
from __future__ import annotations

import argparse, csv, os, random, re, subprocess, sys, time, tempfile
from pathlib import Path

# Stack-limit raise.
try:
    import resource
    _BIG = 1 << 30
    _soft, _hard = resource.getrlimit(resource.RLIMIT_STACK)
    _t = _BIG if _hard == resource.RLIM_INFINITY else min(_BIG, _hard)
    if _t != _soft:
        try: resource.setrlimit(resource.RLIMIT_STACK, (_t, _hard))
        except Exception: resource.setrlimit(resource.RLIMIT_STACK,
                                             (max(_soft, _t-4096), _hard))
except Exception:
    pass

ALGOS = [("Ours_ST", {"PIVOTER_RUN_ST": "1"}), ("REF_R1", {})]


def subsample_edges(in_path: Path, ratio: float, out_path: Path, seed: int = 42):
    """Nested vertex-induced subsampling: reproducibly shuffle the vertex set
    with `seed`, take the first ratio*n vertices, and write the induced
    subgraph (edges with both endpoints kept). Because every call uses the
    same seed and the same input, the 40% sample is a strict superset of
    the 20% sample, etc. Vertex IDs are re-mapped to a contiguous 0..k-1
    range so the header reports |V'| accurately.

    Returns the number of edges kept (for logging)."""
    rng = random.Random(seed)
    with in_path.open() as fin:
        first = fin.readline().split()
        n, m = int(first[0]), int(first[1])
        edges = []
        for line in fin:
            parts = line.split()
            if len(parts) >= 2:
                edges.append((int(parts[0]), int(parts[1])))
    vertices = list(range(n))
    rng.shuffle(vertices)
    keep_n = max(1, int(round(ratio * n)))
    kept = set(vertices[:keep_n])
    # Build contiguous remap so the resulting file has vertex ids in [0, keep_n).
    remap = {v: i for i, v in enumerate(sorted(kept))}
    out_edges = [(remap[u], remap[v]) for (u, v) in edges if u in kept and v in kept]
    with out_path.open("w") as fout:
        fout.write(f"{keep_n} {len(out_edges)}\n")
        for u, v in out_edges:
            fout.write(f"{u} {v}\n")
    return len(out_edges)


def parse_timing(stdout: str):
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', stdout)
    m_mem   = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)', stdout)
    total = float(m_total.group(1)) if m_total else -1.0
    mem   = float(m_mem.group(1))   if m_mem   else -1.0
    return total, mem


def run_one(bin_path: Path, gpath: Path, s: int, env_extra: dict, timeout: int):
    env = os.environ.copy(); env.update(env_extra)
    try:
        proc = subprocess.run([str(bin_path), str(gpath), "1", str(s)],
                              env=env, capture_output=True, text=True,
                              timeout=timeout)
    except subprocess.TimeoutExpired:
        return ("TIMEOUT", -1.0, -1.0)
    if proc.returncode != 0:
        return (f"FAIL({proc.returncode})", -1.0, -1.0)
    t, m = parse_timing(proc.stdout)
    return ("OK" if t >= 0 else "PARSE_FAIL", t, m)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default="./build/bin/degeneracy_cliques")
    ap.add_argument("--graph", required=True, help="path to base .edges")
    ap.add_argument("--ratios", nargs="+", type=float,
                    default=[0.2, 0.4, 0.6, 0.8, 1.0])
    ap.add_argument("--s-list", nargs="+", type=int,
                    default=[2, 3, 5, 7, 9, 11, 13, 15])
    ap.add_argument("--runs", type=int, default=1)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--out", default="paper_data/scalability.csv")
    args = ap.parse_args()

    bin_path = Path(args.bin); gp = Path(args.graph)
    if not bin_path.exists(): sys.exit(f"binary not found: {bin_path}")
    if not gp.exists():       sys.exit(f"graph not found: {gp}")

    out_csv = Path(args.out); out_csv.parent.mkdir(parents=True, exist_ok=True)
    new_file = not out_csv.exists() or out_csv.stat().st_size == 0
    fout = out_csv.open("a")
    if new_file:
        fout.write("graph,base_edges,ratio,kept_edges,s,algorithm,run,time_ms,memory_kB,status\n")
        fout.flush()

    # Parse base edge count.
    with gp.open() as f:
        first = f.readline().split()
        base_edges = int(first[1])
    print(f"[base] {gp.name}  m={base_edges}", flush=True)

    work = Path(tempfile.mkdtemp(prefix="scal_"))
    print(f"[tmp] subsamples in {work}", flush=True)

    try:
        for ratio in args.ratios:
            sub_path = work / f"{gp.stem}_{int(ratio*100):03d}.edges"
            kept = subsample_edges(gp, ratio, sub_path, seed=args.seed)
            print(f"[ratio={ratio:.2f}] kept={kept} edges -> {sub_path.name}",
                  flush=True)

            for s in args.s_list:
                for algo, env_extra in ALGOS:
                    for run in range(args.runs):
                        t0 = time.time()
                        status, t_ms, m_kb = run_one(bin_path, sub_path, s,
                                                    env_extra, args.timeout)
                        elapsed = time.time() - t0
                        t_str = f"{t_ms:.1f}" if t_ms >= 0 else ""
                        m_str = f"{m_kb:.0f}" if m_kb >= 0 else ""
                        fout.write(f"{gp.stem},{base_edges},{ratio},{kept},"
                                   f"{s},{algo},{run},{t_str},{m_str},{status}\n")
                        fout.flush()
                        print(f"  s={s:2d} {algo:7s} run={run} {status:6s} "
                              f"t={t_ms:.0f}ms rss={m_kb:.0f}kB "
                              f"(wall={elapsed:.1f}s)", flush=True)
    finally:
        # Keep tmp around so re-runs don't redo subsampling; comment to clean.
        print(f"[done] subsamples retained in {work}", flush=True)

    fout.close()


if __name__ == "__main__":
    main()
