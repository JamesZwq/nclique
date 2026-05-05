#!/usr/bin/env python3
"""
Stress-test on synthetic dense graphs (parallel, resource-aware).

Reproduces the protocol of Nuclear-CD-TODS §6 (Stress Test), restricted
to r=1.  N=1000, increasing edge density p; for each p run Ours_ST and
REF_R1 across a sweep of s.  Halts when every algorithm has died (i.e.
failed at the lowest s of args.s_list) at some density.

Resource scheduling:
  * MAX_WORKERS hard cap (default 32)
  * Launch gate: do not start a new worker if global used memory exceeds
    MEM_LIMIT_GB
  * Per-process OOM: any single process whose RSS exceeds PER_PROC_MEM_GB
    is killed and recorded as OOM (terminal -- no retry)
  * CPU load-avg gate (--cpu-load-target; disabled by default)
  * Adaptive throttle: longer sleep between launches when memory rises

Usage:
    python3 bench_r1_stress.py \
        --bin ./build/bin/degeneracy_cliques \
        --n 1000 --workers 32 --timeout 1200 \
        --out paper_data/stress.csv
"""
from __future__ import annotations

import argparse, csv, math, multiprocessing as _mp, os, random, re, signal
import subprocess, sys, tempfile, time
from collections import defaultdict
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

# ============ Resource limits ============
# Simple OOM policy: any single process whose RSS exceeds PER_PROC_MEM_GB is
# killed and its (p, s, algo) is recorded as OOM (terminal — no retry).
# A launch gate based on global used memory keeps total demand bounded.
MEM_LIMIT_GB    = 300       # gate: don't launch new worker if used > this
PER_PROC_MEM_GB = 250       # any single proc with RSS > this is OOM-killed
SETTLE_SEC      = 0.2
POLL_SEC        = 3

# Paper §7 calls the static pipeline "Pure" — this is V3 (PIVOTER_RUN_ST_V3),
# the SOTA with fused SDCT + event-driven PeelH.  The previous stress run
# stored Ours_ST (the older PIVOTER_RUN_ST without V2/V3 improvements);
# rerun produces a fresh CSV under a new name so the legacy file stays
# intact for cross-checking.
ALGOS = [("Pure",   {"PIVOTER_RUN_ST_V3": "1"}),
         ("REF_R1", {})]


# ============ Resource probes ============
def get_used_mem_gb() -> float:
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


def get_proc_rss_gb(pid: int) -> float:
    """RSS of a single PID in GB. Returns 0 if process is gone."""
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024 / 1024
    except Exception:
        pass
    return 0.0


def cpu_has_headroom(running_count: int, max_workers: int, cpu_target):
    if running_count >= max_workers:
        return False
    if cpu_target is None:
        return True
    return get_load_avg_1min() < _mp.cpu_count() * cpu_target


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
        return [(i, j) for i in range(n) for j in range(i + 1, n)]
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


# ============ Subprocess output parsing ============
def parse_timing(stdout: str):
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', stdout)
    m_mem   = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)', stdout)
    return (float(m_total.group(1)) if m_total else -1.0,
            float(m_mem.group(1))   if m_mem   else -1.0)


def parse_existing(csv_path: Path):
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
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--densities", nargs="+", type=float,
                    default=[0.01, 0.03, 0.05, 0.08, 0.11, 0.15, 0.20,
                             0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                             0.90, 0.95, 0.99])
    ap.add_argument("--s-list", nargs="+", type=int, default=[5, 7, 9, 11])
    ap.add_argument("--alpha", type=float, default=2.0)
    ap.add_argument("--max-clique", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--timeout", type=int, default=1200)
    ap.add_argument("--cpu-load-target", type=float, default=None)
    ap.add_argument("--out", default="paper_data/15_stress_synthetic_dense_v3.csv",
                    help="default points to the V3-SOTA output CSV; override "
                         "to keep building on the legacy ST CSV.")
    ap.add_argument("--mem-limit-gb", type=float, default=MEM_LIMIT_GB)
    ap.add_argument("--per-proc-mem-gb", type=float, default=PER_PROC_MEM_GB)
    args = ap.parse_args()

    bin_path = Path(args.bin)
    if not bin_path.exists(): sys.exit(f"binary not found: {bin_path}")

    cpu_target = args.cpu_load_target
    max_workers = max(1, min(args.workers, 32))
    mem_limit = args.mem_limit_gb
    per_proc_mem = args.per_proc_mem_gb

    out_csv = Path(args.out); out_csv.parent.mkdir(parents=True, exist_ok=True)
    done = parse_existing(out_csv)
    new_file = not out_csv.exists() or out_csv.stat().st_size == 0

    print(f"[setup] workers={max_workers}  mem_limit={mem_limit:.0f}GB  "
          f"per_proc={per_proc_mem:.0f}GB  timeout={args.timeout}s  "
          f"cpu_target={cpu_target}", flush=True)
    print(f"[setup] {len(done)} (density,s,algo) already OK in {out_csv}",
          flush=True)

    work = Path(tempfile.mkdtemp(prefix="stress_"))
    print(f"[tmp] {work}", flush=True)

    # Generate all graphs upfront.
    print(f"[gen] generating {len(args.densities)} synthetic graphs...",
          flush=True)
    graphs = {}
    for p in args.densities:
        gpath = work / f"stress_n{args.n}_p{int(p*1000):04d}.edges"
        edges = generate_clique_injection(args.n, p, seed=args.seed,
                                          alpha=args.alpha,
                                          max_clique=args.max_clique)
        write_edge_file(gpath, args.n, edges)
        graphs[p] = (gpath, len(edges))
        print(f"  p={p:.3f}  m={len(edges)}  ({gpath.name})", flush=True)

    fout = out_csv.open("a")
    if new_file:
        fout.write("n,density,m,s,algorithm,time_ms,memory_kB,status\n")
        fout.flush()

    # State.
    sfloor = {}             # (p, algo) -> min s that failed
    algo_dead_from = {}     # algo -> min density at which it died
    n_ok = n_failed = n_skipped_floor = n_skipped_dead = n_oom = 0

    def write_row(p, s, algo, status, t, m):
        m_used = graphs[p][1]
        t_str = f"{t:.1f}" if t >= 0 else ""
        m_str = f"{m:.0f}" if m >= 0 else ""
        fout.write(f"{args.n},{p},{m_used},{s},{algo},{t_str},{m_str},{status}\n")
        fout.flush()

    def update_failure(p, s, algo):
        prev = sfloor.get((p, algo), float("inf"))
        if s < prev:
            sfloor[(p, algo)] = s
        if s == args.s_list[0]:
            prev_p = algo_dead_from.get(algo, float("inf"))
            if p < prev_p:
                algo_dead_from[algo] = p
                print(f"[dead] {algo} dead at p>={p} (failed lowest s={s})",
                      flush=True)

    def should_skip(p, s, algo):
        if algo in algo_dead_from and algo_dead_from[algo] <= p:
            return ("SKIP_DEAD_ALGO", -1, -1)
        f = sfloor.get((p, algo))
        if f is not None and s >= f:
            return ("SKIP_FLOOR", -1, -1)
        return None

    def launch(p, s, algo, env_extra) -> dict:
        gpath = graphs[p][0]
        env = os.environ.copy(); env.update(env_extra)
        proc = subprocess.Popen([str(bin_path), str(gpath), "1", str(s)],
                                env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        return dict(proc=proc, p=p, s=s, algo=algo, env=env_extra,
                    t0=time.time(), key=(p, s, algo))

    running = []   # list of dicts above

    def reap():
        nonlocal n_ok, n_failed, n_oom
        finished = []
        for r in running:
            rc = r["proc"].poll()
            if rc is not None:
                finished.append(r)
        for r in finished:
            running.remove(r)
            stdout = r["proc"].stdout.read() if r["proc"].stdout else ""
            r["proc"].stdout.close() if r["proc"].stdout else None
            elapsed = time.time() - r["t0"]
            rc = r["proc"].returncode
            p, s, algo = r["p"], r["s"], r["algo"]
            if rc == 0:
                t_ms, m_kb = parse_timing(stdout)
                if t_ms < 0:
                    write_row(p, s, algo, "PARSE_FAIL", -1, -1)
                    n_failed += 1
                    update_failure(p, s, algo)
                    print(f"  [p={p:.3f} s={s:2d} {algo:7s}]  PARSE_FAIL "
                          f"(wall={elapsed:.1f}s)", flush=True)
                else:
                    write_row(p, s, algo, "OK", t_ms, m_kb)
                    n_ok += 1
                    print(f"  [p={p:.3f} s={s:2d} {algo:7s}]  OK       "
                          f"t={t_ms:.0f}ms rss={m_kb:.0f}kB "
                          f"(wall={elapsed:.1f}s)", flush=True)
            else:
                # rc != 0: could be OOM kill (-9 / -SIGKILL), other fail, or
                # timeout (we don't enforce timeout here yet).
                status = f"FAIL({rc})"
                write_row(p, s, algo, status, -1, -1)
                n_failed += 1
                update_failure(p, s, algo)
                print(f"  [p={p:.3f} s={s:2d} {algo:7s}]  {status:12s} "
                      f"(wall={elapsed:.1f}s)", flush=True)

    def check_timeouts():
        now = time.time()
        to_kill = [r for r in running if (now - r["t0"]) > args.timeout]
        for r in to_kill:
            try:
                r["proc"].kill()
                r["proc"].wait(timeout=10)
            except Exception: pass
            running.remove(r)
            p, s, algo = r["p"], r["s"], r["algo"]
            write_row(p, s, algo, "TIMEOUT", -1, -1)
            update_failure(p, s, algo)
            print(f"  [p={p:.3f} s={s:2d} {algo:7s}]  TIMEOUT "
                  f"(wall={time.time()-r['t0']:.0f}s)", flush=True)

    def check_proc_mem():
        """Per-process OOM gate: any single proc with RSS > per_proc_mem is
        terminally killed and its (p, s, algo) recorded as OOM (no retry)."""
        nonlocal n_oom
        for r in list(running):
            rss = get_proc_rss_gb(r["proc"].pid)
            if rss > per_proc_mem:
                try:
                    r["proc"].kill(); r["proc"].wait(timeout=10)
                except Exception: pass
                running.remove(r)
                p, s, algo = r["p"], r["s"], r["algo"]
                write_row(p, s, algo, "OOM", -1, -1)
                n_oom += 1
                update_failure(p, s, algo)
                print(f"  [p={p:.3f} s={s:2d} {algo:7s}]  OOM "
                      f"RSS={rss:.0f}GB > {per_proc_mem:.0f}GB", flush=True)

    # Build job list.
    all_jobs = []
    for p in args.densities:
        for s in args.s_list:
            for algo, env in ALGOS:
                if (f"{p}", s, algo) in done:
                    continue
                all_jobs.append((p, s, algo, env))

    print(f"[setup] {len(all_jobs)} jobs to attempt (after resume)", flush=True)

    # ===== Main loop =====
    ji = 0
    t_start = time.time()
    halt_announced = False
    while ji < len(all_jobs) or running:
        reap()
        check_timeouts()
        check_proc_mem()

        # If every algo is dead, drain and stop.
        if (not halt_announced
            and all(a in algo_dead_from for a, _ in ALGOS)):
            print(f"[stop] all algorithms dead; draining {len(running)} "
                  f"in-flight jobs", flush=True)
            halt_announced = True
            ji = len(all_jobs)

        # Pick next job from the main queue.
        job = None
        if ji < len(all_jobs):
            job = all_jobs[ji]
            ji += 1
        else:
            if running: time.sleep(POLL_SEC)
            continue

        p, s, algo, env_extra = job

        # Skip checks.
        skip = should_skip(p, s, algo)
        if skip is not None:
            status, _, _ = skip
            write_row(p, s, algo, status, -1, -1)
            if status == "SKIP_DEAD_ALGO": n_skipped_dead += 1
            else: n_skipped_floor += 1
            continue

        # Wait for capacity.
        wait_start = time.time()
        while ((not cpu_has_headroom(len(running), max_workers, cpu_target))
               or get_used_mem_gb() >= mem_limit):
            reap(); check_timeouts(); check_proc_mem()
            time.sleep(POLL_SEC)
            # Status print every minute.
            if int(time.time() - wait_start) % 60 == 0 and len(running) > 0:
                mem = get_used_mem_gb()
                load = get_load_avg_1min()
                print(f"[wait] running={len(running)} pending={len(all_jobs)-ji} "
                      f"mem={mem:.0f}GB load={load:.1f}", flush=True)

        # Re-check skip after waiting.
        skip = should_skip(p, s, algo)
        if skip is not None:
            status, _, _ = skip
            write_row(p, s, algo, status, -1, -1)
            if status == "SKIP_DEAD_ALGO": n_skipped_dead += 1
            else: n_skipped_floor += 1
            continue

        # Launch.
        running.append(launch(p, s, algo, env_extra))

        # Adaptive throttle: longer pauses when memory is climbing.
        mem_now = get_used_mem_gb()
        if   mem_now > mem_limit * 0.85: time.sleep(3)
        elif mem_now > mem_limit * 0.5:  time.sleep(1)
        else:                             time.sleep(SETTLE_SEC)

        # Periodic scheduler status.
        if (ji % 10 == 0 or len(running) >= max_workers):
            mem = get_used_mem_gb(); load = get_load_avg_1min()
            print(f"[sched] running={len(running):2d}  pending={len(all_jobs)-ji:3d} "
                  f"ok={n_ok} fail={n_failed} oom={n_oom} "
                  f"skip(floor/dead)={n_skipped_floor}/{n_skipped_dead}  "
                  f"mem={mem:.0f}GB load={load:.1f}  "
                  f"elapsed={time.time()-t_start:.0f}s", flush=True)

    # Drain remaining.
    print(f"[drain] {len(running)} jobs still running", flush=True)
    while running:
        reap(); check_timeouts(); check_proc_mem()
        time.sleep(POLL_SEC)

    fout.close()
    print(f"\n=== final summary ===", flush=True)
    print(f"  ok={n_ok}  failed={n_failed}  oom={n_oom}  "
          f"skip(floor/dead)={n_skipped_floor}/{n_skipped_dead}", flush=True)
    print(f"  algo_dead_from = {algo_dead_from}", flush=True)
    print(f"  elapsed = {time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
