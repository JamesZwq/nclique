#!/usr/bin/env python3
"""
Stress-test on synthetic dense graphs: incremental clique-injection.

Reproduces the protocol of Nuclear-CD-TODS §6 (Stress Test), restricted
to r=1.  N=1000, increasing edge density p; for each p run Ours_ST and
REF_R1 across a sweep of s.

Graph generator: incremental clique injection.
    - sample a clique size from a truncated power-law (max 100)
    - sample that many distinct vertices uniformly
    - add the all-pairs edges of the clique
    - repeat until the desired edge count is reached
This produces graphs with the heavy clique structure that maximises
the SOTA's combinatorial blow-up.

Usage:
    python3 bench_r1_stress.py \
        --bin ./build/bin/degeneracy_cliques \
        --n 1000 --densities 0.01 0.03 0.05 0.08 0.11 0.15 0.20 \
        --s-list 5 7 9 11 \
        --out paper_data/stress.csv
"""
from __future__ import annotations

import argparse, csv, math, os, random, re, subprocess, sys, tempfile, time
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


def truncated_powerlaw(rng: random.Random, alpha: float, lo: int, hi: int) -> int:
    """Sample a positive integer from a truncated power law on [lo, hi]
    with exponent alpha (P(k) ~ k^-alpha).  Inverse-transform sampling."""
    u = rng.random()
    # CDF F(k) = (k^(1-alpha) - lo^(1-alpha)) / (hi^(1-alpha) - lo^(1-alpha))
    # Solve F(k) = u for k.
    a = 1.0 - alpha
    lo_a = lo ** a; hi_a = hi ** a
    val = (lo_a + u * (hi_a - lo_a)) ** (1.0 / a)
    k = int(round(val))
    return max(lo, min(hi, k))


def generate_clique_injection(n: int, p: float, seed: int = 42,
                              alpha: float = 2.0, max_clique: int = 100) -> list:
    """Generate a graph by injecting cliques until edge density >= p.
    Returns list of (u, v) edges (u < v)."""
    rng = random.Random(seed)
    target_edges = int(round(p * n * (n - 1) / 2))
    edges = set()
    while len(edges) < target_edges:
        k = truncated_powerlaw(rng, alpha, 3, min(max_clique, n))
        nodes = rng.sample(range(n), k)
        for i in range(k):
            for j in range(i + 1, k):
                u, v = nodes[i], nodes[j]
                if u > v: u, v = v, u
                edges.add((u, v))
                if len(edges) >= target_edges: break
            if len(edges) >= target_edges: break
    return list(edges)


def write_edge_file(path: Path, n: int, edges: list):
    with path.open("w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")


def parse_timing(stdout: str):
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', stdout)
    m_mem   = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)', stdout)
    return (float(m_total.group(1)) if m_total else -1.0,
            float(m_mem.group(1))   if m_mem   else -1.0)


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
    ap.add_argument("--n", type=int, default=1000,
                    help="vertex count of synthetic graphs")
    ap.add_argument("--densities", nargs="+", type=float,
                    default=[0.01, 0.03, 0.05, 0.08, 0.11, 0.15, 0.20])
    ap.add_argument("--s-list", nargs="+", type=int,
                    default=[5, 7, 9, 11])
    ap.add_argument("--alpha", type=float, default=2.0,
                    help="power-law exponent for clique-size sampling")
    ap.add_argument("--max-clique", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--out", default="paper_data/stress.csv")
    args = ap.parse_args()

    bin_path = Path(args.bin)
    if not bin_path.exists(): sys.exit(f"binary not found: {bin_path}")

    out_csv = Path(args.out); out_csv.parent.mkdir(parents=True, exist_ok=True)
    new_file = not out_csv.exists() or out_csv.stat().st_size == 0
    fout = out_csv.open("a")
    if new_file:
        fout.write("n,density,m,s,algorithm,time_ms,memory_kB,status\n")
        fout.flush()

    work = Path(tempfile.mkdtemp(prefix="stress_"))
    print(f"[tmp] {work}", flush=True)

    timeout_floor = {}   # (density, algo) -> int s
    for p in args.densities:
        gpath = work / f"stress_n{args.n}_p{int(p*1000):04d}.edges"
        edges = generate_clique_injection(args.n, p, seed=args.seed,
                                          alpha=args.alpha,
                                          max_clique=args.max_clique)
        write_edge_file(gpath, args.n, edges)
        m = len(edges)
        print(f"[p={p:.2f}] n={args.n} m={m} (target {int(p*args.n*(args.n-1)/2)})",
              flush=True)

        for s in args.s_list:
            for algo, env_extra in ALGOS:
                key = (p, algo)
                floor = timeout_floor.get(key)
                if floor is not None and s >= floor:
                    fout.write(f"{args.n},{p},{m},{s},{algo},,,SKIP_FLOOR\n")
                    fout.flush()
                    print(f"  s={s:2d} {algo:7s}  SKIP (floor at s={floor})",
                          flush=True)
                    continue
                t0 = time.time()
                status, t_ms, m_kb = run_one(bin_path, gpath, s, env_extra,
                                              args.timeout)
                elapsed = time.time() - t0
                t_str = f"{t_ms:.1f}" if t_ms >= 0 else ""
                m_str = f"{m_kb:.0f}" if m_kb >= 0 else ""
                fout.write(f"{args.n},{p},{m},{s},{algo},{t_str},{m_str},{status}\n")
                fout.flush()
                print(f"  s={s:2d} {algo:7s}  {status:8s} "
                      f"t={t_ms:.0f}ms rss={m_kb:.0f}kB (wall={elapsed:.1f}s)",
                      flush=True)
                if status in ("TIMEOUT",) or status.startswith("FAIL"):
                    timeout_floor[key] = s
                    print(f"  -> floor: ({p}, {algo}) s>={s}", flush=True)

    fout.close()
    print(f"[done] {out_csv}", flush=True)


if __name__ == "__main__":
    main()
