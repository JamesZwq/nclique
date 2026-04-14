#!/usr/bin/env python3
"""
V3 Benchmark: ST vs V3 vs V3_NP across all (r,s) combinations.
Auto-probes max clique size. Skips if algo already timed out for smaller r at same s.

Usage:
  nohup python3 bench_v3_all.py > bench_v3_all.log 2>&1 &
"""

import subprocess, os, sys, time, re, csv, json
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

# ============ Config ============
BIN = "./build/bin/degeneracy_cliques"
TIMEOUT = 600  # seconds per job
MAX_WORKERS = 8
OUTCSV = "bench_v3_all_results.csv"
LOGDIR = Path("bench_v3_all_logs")
DATADIR = "/data/wenqianz"

GRAPHS = ["com-dblp", "web-Stanford", "dblp-core30", "email-Eu-core",
          "com-youtube", "web-it-2004"]

ALGOS = {
    "ST":    {"env": "PIVOTER_RUN_ST"},
    "V3":    {"env": "PIVOTER_RUN_REGION_V3"},
    "V3_NP": {"env": "PIVOTER_RUN_REGION_V3", "extra": {"PIVOTER_V3_NO_PRIVATE": "1"}},
}

# ============ Helpers ============
def link_graphs():
    os.makedirs("graphs", exist_ok=True)
    for g in GRAPHS:
        f = f"graphs/{g}.edges"
        src = f"{DATADIR}/{g}.edges"
        if not os.path.exists(f) and os.path.exists(src):
            os.symlink(src, f)
            print(f"  Linked {g}.edges")

def build():
    print("Building...")
    os.makedirs("build", exist_ok=True)
    subprocess.run("cmake -S . -B build -DCMAKE_BUILD_TYPE=Release".split(),
                   capture_output=True)
    r = subprocess.run("cmake --build build -j 12 --target degeneracy_cliques".split(),
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"Build failed:\n{r.stderr[-500:]}")
        sys.exit(1)
    print("  Build OK")

def probe_max_clique(graph):
    """Find max clique size by running with decreasing s until MCs found."""
    gf = f"graphs/{graph}.edges"
    if not os.path.exists(gf):
        return 0

    # First get a rough upper bound from MaxCliqEnum with s=4
    try:
        out = subprocess.run(
            [BIN, gf, "3", "4"],
            capture_output=True, text=True, timeout=120,
            env={**os.environ, "PIVOTER_RUN_REGION_V2": "1"}
        )
        txt = out.stdout + out.stderr
        # Look for "minSize=K" in MaxCliqEnum line
        m = re.search(r'minSize=(\d+)', txt)
        if m:
            min_size = int(m.group(1))
        else:
            return 4  # fallback
    except (subprocess.TimeoutExpired, Exception):
        return 4

    # Binary search for max s where MCs exist
    lo, hi = min_size, 300
    while lo < hi:
        mid = (lo + hi + 1) // 2
        try:
            out = subprocess.run(
                [BIN, gf, "3", str(mid)],
                capture_output=True, text=True, timeout=30,
                env={**os.environ, "PIVOTER_RUN_REGION_V2": "1"}
            )
            txt = out.stdout + out.stderr
            m = re.search(r'Maximal cliques: (\d+)', txt)
            if m and int(m.group(1)) > 0:
                lo = mid
            else:
                hi = mid - 1
        except (subprocess.TimeoutExpired, Exception):
            hi = mid - 1

    return lo

def run_job(graph, r, s, algo_name, algo_cfg):
    """Run a single benchmark job. Returns (graph, r, s, algo, status, log_text)."""
    gf = f"graphs/{graph}.edges"
    logfile = LOGDIR / f"{graph}_r{r}_s{s}_{algo_name}.log"

    env = {**os.environ, algo_cfg["env"]: "1"}
    if "extra" in algo_cfg:
        env.update(algo_cfg["extra"])

    try:
        result = subprocess.run(
            [BIN, gf, str(r), str(s)],
            capture_output=True, text=True, timeout=TIMEOUT, env=env
        )
        txt = result.stdout + result.stderr
        status = "OK" if result.returncode == 0 else f"ERROR({result.returncode})"
    except subprocess.TimeoutExpired:
        txt = "TIMEOUT"
        status = "TIMEOUT"

    logfile.write_text(txt)

    # Extract timing
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', txt)
    m_peel = re.search(r'Peeling time:\s*([\d.]+)', txt)
    total_ms = float(m_total.group(1)) if m_total else -1
    peel_ms = float(m_peel.group(1)) if m_peel else -1

    return graph, r, s, algo_name, status, total_ms, peel_ms

def load_existing():
    """Load already-completed results to skip."""
    done = set()
    if os.path.exists(OUTCSV):
        with open(OUTCSV) as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = (row["graph"], int(row["r"]), int(row["s"]), row["algo"])
                done.add(key)
    return done

# ============ Main ============
def main():
    print("=" * 60)
    print("  V3 Full Benchmark")
    print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    link_graphs()
    build()
    LOGDIR.mkdir(exist_ok=True)

    # Probe max clique sizes
    print("\nProbing max clique sizes...")
    max_cliques = {}
    for graph in GRAPHS:
        if not os.path.exists(f"graphs/{graph}.edges"):
            print(f"  {graph}: SKIP (not found)")
            continue
        mc = probe_max_clique(graph)
        max_cliques[graph] = mc
        print(f"  {graph}: max_clique_size={mc}")

    # Generate all jobs, ordered by (graph, s, r) for timeout skipping
    all_jobs = []
    for graph in GRAPHS:
        if graph not in max_cliques:
            continue
        max_k = max_cliques[graph]
        for s in range(4, max_k + 1):
            for r in range(3, s):
                for algo_name in ALGOS:
                    all_jobs.append((graph, r, s, algo_name))

    done = load_existing()
    print(f"\nTotal (r,s,algo) combinations: {len(all_jobs)}")
    print(f"Already completed: {len(done)}")

    # Setup CSV
    fieldnames = ["graph", "r", "s", "algo", "status", "total_ms", "peel_ms"]
    if not os.path.exists(OUTCSV):
        with open(OUTCSV, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=fieldnames).writeheader()

    # Timeout tracking: per (graph, algo, s) → min r that timed out
    # If algo X timed out at (r, s), skip (r', s) for r' > r
    timeout_at = defaultdict(lambda: float('inf'))  # (graph, algo, s) → min_timeout_r

    # Process jobs in order (by graph, then s ascending, then r ascending)
    # Run parallel within each (graph, s) group
    launched = 0
    skipped_timeout = 0

    # Process one graph at a time. Within each graph, process s ascending.
    # For each s, run all (r, algo) in parallel but check timeout skips.
    # Strategy: for each s, launch all valid (r, algo) jobs in parallel.
    # After batch completes, update timeout_at for next s.
    for graph in GRAPHS:
        if graph not in max_cliques:
            continue
        max_k = max_cliques[graph]
        print(f"\n--- {graph} (max_clique={max_k}) ---")

        for s in range(4, max_k + 1):
            # Collect jobs for this (graph, s)
            batch = []
            for r in range(3, s):
                for algo_name in ALGOS:
                    key = (graph, r, s, algo_name)
                    if key in done:
                        continue
                    min_to_r = timeout_at[(graph, algo_name, s)]
                    if r >= min_to_r:
                        skipped_timeout += 1
                        with open(OUTCSV, "a", newline="") as f:
                            w = csv.DictWriter(f, fieldnames=fieldnames)
                            w.writerow({"graph": graph, "r": r, "s": s,
                                        "algo": algo_name, "status": "SKIP_TIMEOUT",
                                        "total_ms": "", "peel_ms": ""})
                        done.add(key)
                        continue
                    batch.append((graph, r, s, algo_name))

            if not batch:
                continue

            # Run batch in parallel
            with ProcessPoolExecutor(max_workers=min(MAX_WORKERS, len(batch))) as ex:
                futures = {}
                for g, r, s2, an in batch:
                    f = ex.submit(run_job, g, r, s2, an, ALGOS[an])
                    futures[f] = (g, r, s2, an)

                for f in as_completed(futures):
                    g, rr, ss, an, status, total_ms, peel_ms = f.result()
                    launched += 1
                    t_str = f"{total_ms:.0f}ms" if total_ms > 0 else "N/A"
                    print(f"  {an:>6} {g} r={rr} s={ss} {status} total={t_str}")

                    with open(OUTCSV, "a", newline="") as f2:
                        w = csv.DictWriter(f2, fieldnames=fieldnames)
                        w.writerow({"graph": g, "r": rr, "s": ss, "algo": an,
                                    "status": status,
                                    "total_ms": f"{total_ms:.1f}" if total_ms > 0 else "",
                                    "peel_ms": f"{peel_ms:.1f}" if peel_ms > 0 else ""})
                    done.add((g, rr, ss, an))

                    if status in ("TIMEOUT", "OOM"):
                        timeout_at[(g, an, ss)] = min(timeout_at[(g, an, ss)], rr)
                        # Also propagate: if timed out at (r, s), skip (r, s') for s' > s
                        # (larger s with same r is likely harder)
                        for s_future in range(ss + 1, max_cliques.get(g, 0) + 1):
                            timeout_at[(g, an, s_future)] = min(
                                timeout_at[(g, an, s_future)], rr)

    print(f"\nDone. Launched: {launched}, skipped (timeout): {skipped_timeout}")
    print(f"Results: {OUTCSV}")
    print(f"Logs: {LOGDIR}/")

if __name__ == "__main__":
    main()
