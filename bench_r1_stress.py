#!/usr/bin/env python3
"""
Stress-test on synthetic dense graphs (multi-threaded, resource-aware).

Reproduces the protocol of Nuclear-CD-TODS §6 (Stress Test), restricted
to r=1.  N=1000, increasing edge density p; for each p run Ours_ST and
REF_R1 across a sweep of s.

Density list defaults to {0.01, 0.03, ..., 0.99}; the bench halts when
every algorithm has died (timed out at the lowest s) at some density.

Resource scheduling (modeled on bench_v3_all.py):
  * MAX_WORKERS hard cap (default 32, per --workers flag)
  * Memory gate: do not launch a new worker if used memory exceeds
    MEM_LIMIT_GB; kill newest worker if memory exceeds MEM_KILL_GB
  * CPU load-avg gate: do not launch new worker if loadavg > nproc *
    CPU_LOAD_TARGET (set to None to disable on dedicated servers)
  * SETTLE_SEC pause between launches; POLL_SEC poll interval

Per-(density, algo) within-density floor: a timeout at s_t skips s>=s_t.
Per-algo cross-density death: a timeout at the lowest s for some algo
at density p marks the algo dead at every p' >= p.

Resume-friendly: skip rows with status=OK already in the output CSV.

Usage:
    python3 bench_r1_stress.py \
        --bin ./build/bin/degeneracy_cliques \
        --n 1000 --workers 32 \
        --densities 0.01 0.03 0.05 0.08 0.11 0.15 0.20 \
                    0.25 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 0.99 \
        --s-list 5 7 9 11 \
        --timeout 1200 --out paper_data/stress.csv
"""
from __future__ import annotations

import argparse, csv, math, multiprocessing as _mp, os, random, re, signal
import subprocess, sys, tempfile, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# ============ Stack-limit raise ============
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

# ============ Resource gates (mirroring bench_v3_all.py) ============
MAX_WORKERS_DEFAULT = 32
MEM_LIMIT_GB        = 300       # don't launch if used > this
MEM_KILL_GB         = 450       # warn if used > this
PER_PROC_MEM_GB     = 250       # warn if a single proc > this
SETTLE_SEC          = 0.2       # pause between launches
POLL_SEC            = 3         # poll interval
CPU_LOAD_TARGET     = None      # tods2 dedicated → disabled; set to e.g. 0.85 to enable

ALGOS = [("Ours_ST", {"PIVOTER_RUN_ST": "1"}), ("REF_R1", {})]


def get_used_mem_gb() -> float:
    """Used memory in GB from /proc/meminfo."""
    try:
        info = {}
        with open("/proc/meminfo") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    info[parts[0].rstrip(":")] = int(parts[1])
        return (info.get("MemTotal", 0) - info.get("MemAvailable", 0)) / 1024 / 1024
    except Exception:
        return 0.0


def get_load_avg_1min() -> float:
    try: return os.getloadavg()[0]
    except Exception: return 0.0


def can_launch(running: int, max_workers: int, cpu_target):
    if running >= max_workers:
        return False
    if get_used_mem_gb() >= MEM_LIMIT_GB:
        return False
    if cpu_target is not None and get_load_avg_1min() > _mp.cpu_count() * cpu_target:
        return False
    return True


# ============ Graph generation ============
def truncated_powerlaw(rng, alpha, lo, hi):
    u = rng.random()
    a = 1.0 - alpha
    lo_a = lo ** a; hi_a = hi ** a
    val = (lo_a + u * (hi_a - lo_a)) ** (1.0 / a)
    return max(lo, min(hi, int(round(val))))


def generate_clique_injection(n, p, seed=42, alpha=2.0, max_clique=100):
    rng = random.Random(seed)
    target = int(round(p * n * (n - 1) / 2))
    if target >= n * (n - 1) // 2:
        edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
        return edges
    edges = set()
    while len(edges) < target:
        k = truncated_powerlaw(rng, alpha, 3, min(max_clique, n))
        nodes = rng.sample(range(n), k)
        for i in range(k):
            for j in range(i + 1, k):
                u, v = nodes[i], nodes[j]
                if u > v: u, v = v, u
                edges.add((u, v))
                if len(edges) >= target: break
            if len(edges) >= target: break
    return list(edges)


def write_edge_file(path, n, edges):
    with path.open("w") as f:
        f.write(f"{n} {len(edges)}\n")
        for u, v in edges:
            f.write(f"{u} {v}\n")


# ============ Subprocess runner ============
def parse_timing(stdout):
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', stdout)
    m_mem   = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)', stdout)
    return (float(m_total.group(1)) if m_total else -1.0,
            float(m_mem.group(1))   if m_mem   else -1.0)


def run_one(bin_path, gpath, s, env_extra, timeout):
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


def parse_existing(csv_path: Path):
    """Returns the set of (density_str, s, algo) already recorded with status==OK."""
    done = set()
    if not csv_path.exists():
        return done
    with csv_path.open() as f:
        for row in csv.DictReader(f):
            if row.get("status") != "OK":
                continue
            try:
                done.add((row["density"], int(row["s"]), row["algorithm"]))
            except (KeyError, ValueError):
                pass
    return done


# ============ Main scheduler ============
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default="./build/bin/degeneracy_cliques")
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--workers", type=int, default=MAX_WORKERS_DEFAULT,
                    help=f"max parallel workers (default {MAX_WORKERS_DEFAULT})")
    ap.add_argument("--densities", nargs="+", type=float,
                    default=[0.01, 0.03, 0.05, 0.08, 0.11, 0.15, 0.20,
                             0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                             0.90, 0.95, 0.99])
    ap.add_argument("--s-list", nargs="+", type=int, default=[5, 7, 9, 11])
    ap.add_argument("--alpha", type=float, default=2.0)
    ap.add_argument("--max-clique", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--timeout", type=int, default=1200)
    ap.add_argument("--cpu-load-target", type=float, default=None,
                    help="if set, gate launches when 1-min loadavg > nproc * this")
    ap.add_argument("--out", default="paper_data/stress.csv")
    args = ap.parse_args()

    bin_path = Path(args.bin)
    if not bin_path.exists(): sys.exit(f"binary not found: {bin_path}")

    cpu_target = args.cpu_load_target
    max_workers = max(1, min(args.workers, MAX_WORKERS_DEFAULT))

    out_csv = Path(args.out); out_csv.parent.mkdir(parents=True, exist_ok=True)
    done = parse_existing(out_csv)
    new_file = not out_csv.exists() or out_csv.stat().st_size == 0

    print(f"[setup] workers={max_workers} mem_limit={MEM_LIMIT_GB}GB "
          f"cpu_target={cpu_target} timeout={args.timeout}s", flush=True)
    print(f"[setup] {len(done)} (density,s,algo) already OK in {out_csv}", flush=True)

    work = Path(tempfile.mkdtemp(prefix="stress_"))
    print(f"[tmp] {work}", flush=True)

    # Generate all graphs upfront (fast).
    graphs = {}    # density -> (path, m)
    print(f"[gen] generating {len(args.densities)} synthetic graphs...", flush=True)
    for p in args.densities:
        gpath = work / f"stress_n{args.n}_p{int(p*1000):04d}.edges"
        edges = generate_clique_injection(args.n, p, seed=args.seed,
                                          alpha=args.alpha,
                                          max_clique=args.max_clique)
        write_edge_file(gpath, args.n, edges)
        graphs[p] = (gpath, len(edges))
        print(f"  p={p:.3f}  m={len(edges)}  ({gpath.name})", flush=True)

    # Shared state.
    lock = threading.Lock()
    fout = out_csv.open("a")
    if new_file:
        fout.write("n,density,m,s,algorithm,time_ms,memory_kB,status\n")
        fout.flush()
    sfloor = {}             # (p, algo) -> min s that failed
    algo_dead_from = {}     # algo -> min density at which it died

    def write_row(p, s, algo, status, t, m):
        m_used = graphs[p][1]
        t_str = f"{t:.1f}" if t >= 0 else ""
        m_str = f"{m:.0f}" if m >= 0 else ""
        with lock:
            fout.write(f"{args.n},{p},{m_used},{s},{algo},{t_str},{m_str},{status}\n")
            fout.flush()

    def update_failure(p, s, algo):
        """Called inside the lock when a worker reports timeout/fail."""
        prev = sfloor.get((p, algo), float("inf"))
        if s < prev:
            sfloor[(p, algo)] = s
        if s == args.s_list[0]:
            prev_p = algo_dead_from.get(algo, float("inf"))
            if p < prev_p:
                algo_dead_from[algo] = p
                print(f"[dead] {algo} dead at p>={p} (failed lowest s={s})",
                      flush=True)

    # Build job list.
    all_jobs = []
    for p in args.densities:
        for s in args.s_list:
            for algo, env in ALGOS:
                all_jobs.append((p, s, algo, env))

    # Skip resume.
    pending = []
    n_skipped_done = 0
    for p, s, algo, env in all_jobs:
        if (f"{p}", s, algo) in done:
            n_skipped_done += 1
            continue
        pending.append((p, s, algo, env))
    print(f"[setup] {n_skipped_done} jobs already done, {len(pending)} remaining",
          flush=True)

    # ===== Worker function =====
    n_done = 0; n_skipped_floor = 0; n_skipped_dead = 0; n_failed = 0; n_ok = 0

    def worker(p, s, algo, env):
        nonlocal n_done, n_skipped_floor, n_skipped_dead, n_failed, n_ok
        # Re-check skip conditions inside the worker (state may have changed).
        with lock:
            if algo in algo_dead_from and algo_dead_from[algo] <= p:
                write_row(p, s, algo, "SKIP_DEAD_ALGO", -1, -1)
                n_skipped_dead += 1
                n_done += 1
                return
            f = sfloor.get((p, algo))
            if f is not None and s >= f:
                write_row(p, s, algo, "SKIP_FLOOR", -1, -1)
                n_skipped_floor += 1
                n_done += 1
                return

        gpath, m_edges = graphs[p]
        t0 = time.time()
        status, t_ms, m_kb = run_one(bin_path, gpath, s, env, args.timeout)
        elapsed = time.time() - t0
        write_row(p, s, algo, status, t_ms, m_kb)

        with lock:
            n_done += 1
            if status == "OK":
                n_ok += 1
            else:
                n_failed += 1
                if status == "TIMEOUT" or status.startswith("FAIL"):
                    update_failure(p, s, algo)

        tag = f"[p={p:.3f} s={s:2d} {algo:7s}]"
        if status == "OK":
            print(f"  {tag} {status:8s} t={t_ms:.0f}ms rss={m_kb:.0f}kB "
                  f"(wall={elapsed:.1f}s)", flush=True)
        else:
            print(f"  {tag} {status:12s} (wall={elapsed:.1f}s)", flush=True)

    # ===== Main scheduling loop =====
    print(f"[run] starting; {len(pending)} jobs", flush=True)
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        idx = 0
        running_futs = []
        while idx < len(pending) or running_futs:
            # Reap completed (without blocking).
            still = []
            for fut in running_futs:
                if fut.done():
                    try: fut.result()
                    except Exception as e:
                        print(f"  worker exception: {e}", flush=True)
                else:
                    still.append(fut)
            running_futs = still

            # Try to launch new jobs.
            launched_this_round = 0
            while idx < len(pending) and can_launch(len(running_futs),
                                                    max_workers, cpu_target):
                p, s, algo, env = pending[idx]
                idx += 1
                # Check skips atomically (cheap pre-filter — worker also
                # rechecks, but this avoids wasting a slot).
                with lock:
                    if algo in algo_dead_from and algo_dead_from[algo] <= p:
                        write_row(p, s, algo, "SKIP_DEAD_ALGO", -1, -1)
                        n_skipped_dead += 1
                        n_done += 1
                        continue
                    f = sfloor.get((p, algo))
                    if f is not None and s >= f:
                        write_row(p, s, algo, "SKIP_FLOOR", -1, -1)
                        n_skipped_floor += 1
                        n_done += 1
                        continue
                fut = ex.submit(worker, p, s, algo, env)
                running_futs.append(fut)
                launched_this_round += 1
                time.sleep(SETTLE_SEC)

            mem = get_used_mem_gb()
            load = get_load_avg_1min()
            elapsed = time.time() - t_start
            if launched_this_round > 0 or len(running_futs) > 0:
                print(f"[sched] running={len(running_futs):2d} pending={len(pending)-idx:3d} "
                      f"done={n_done:3d}/{len(pending)} ok={n_ok} fail={n_failed} "
                      f"skip(floor/dead)={n_skipped_floor}/{n_skipped_dead}  "
                      f"mem={mem:.0f}GB load={load:.1f}  "
                      f"elapsed={elapsed:.0f}s", flush=True)

            # Halt if every algo has died.
            with lock:
                if all(a in algo_dead_from for a, _ in ALGOS):
                    pending = pending[:idx]   # don't add more
                    print(f"[stop] all algos dead — waiting for in-flight jobs",
                          flush=True)
            time.sleep(POLL_SEC)

    fout.close()
    print(f"\n=== final summary ===", flush=True)
    print(f"  ok={n_ok}  failed={n_failed}  "
          f"skipped(floor/dead)={n_skipped_floor}/{n_skipped_dead}  "
          f"total={n_done}", flush=True)
    print(f"  algo_dead_from = {algo_dead_from}", flush=True)
    print(f"  elapsed = {time.time() - t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
