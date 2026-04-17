#!/usr/bin/env python3
"""
V3 Benchmark: ST vs V3 vs V3_NP across all (r,s) combinations.
Auto-probes max clique size. Memory-aware parallel scheduling.
Timeout skip: if algo times out at (r,s), skip r'>=r for s'>=s.

Usage:
  nohup python3 bench_v3_all.py > bench_v3_all.log 2>&1 &
"""

import subprocess, os, sys, time, re, csv, signal, json
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

# Raise stack limit so deep BK recursion (e.g., high-degree hubs in soc-Epinions1)
# doesn't segfault. Linux default is 8MB; we raise as high as the hard limit allows.
try:
    import resource
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    # macOS rejects if target == hard exactly; use hard-1 as a safe upper bound.
    if hard == resource.RLIM_INFINITY:
        target = resource.RLIM_INFINITY
    else:
        target = max(soft, hard - 4096)
    if target > soft:
        resource.setrlimit(resource.RLIMIT_STACK, (target, hard))
    new_soft = resource.getrlimit(resource.RLIMIT_STACK)[0]
    def _fmt(x):
        return "unlimited" if x == resource.RLIM_INFINITY else f"{x/1024/1024:.1f}MB"
    print(f"[stack] RLIMIT_STACK: {_fmt(soft)} -> {_fmt(new_soft)} (hard={_fmt(hard)})", flush=True)
except Exception as e:
    print(f"[stack] WARNING: failed to raise stack limit: {e}", flush=True)

# ============ Config ============
BIN = "./build/bin/degeneracy_cliques"
TIMEOUT = 3600         # seconds per job (1 hour)
MAX_WORKERS = 32
MEM_LIMIT_GB = 300     # don't launch if total used > this
MEM_KILL_GB = 450      # kill newest if total used > this
PER_PROC_MEM_GB = 250  # kill individual process if RSS > this
SETTLE_SEC = 0.1       # brief pause between launches
POLL_SEC = 3           # poll interval
OUTCSV = "bench_v3_all_results.csv"
LOGDIR = Path("bench_v3_all_logs")
DATADIR = "/data/wenqianz"
# Split graphs across servers: tods1 gets dense graphs, tods2 gets large/power-law
SERVER_GRAPHS = {
    "tods1": ["com-dblp", "email-Eu-core", "com-youtube", "soc-Epinions1"],
    "tods2": ["dblp-core30", "web-Stanford", "web-it-2004"],
}



ALL_GRAPHS = SERVER_GRAPHS["tods1"] + SERVER_GRAPHS["tods2"]


import socket
def get_graphs():
    """Pick graphs based on hostname or CLI arg."""
    # CLI: python3 bench_v3_all.py tods1
    if len(sys.argv) > 1 and sys.argv[1] in SERVER_GRAPHS:
        return SERVER_GRAPHS[sys.argv[1]]
    # Auto-detect from hostname
    hostname = socket.gethostname().lower()
    for key in SERVER_GRAPHS:
        if key in hostname:
            return SERVER_GRAPHS[key]
    # Default: all graphs
    return ALL_GRAPHS

GRAPHS = get_graphs()

ALGOS = {
    "REF":       {"env": "PIVOTER_RUN_REF"},
    "ST":        {"env": "PIVOTER_RUN_ST"},
    # "V3Fast" = current optimized V3 with Private Cloud ON. The legacy "V3"
    # and "V3_NP" keys (slow, sometimes-wrong due to uint32 overflow) have
    # been retired; older rows under those keys remain in the CSV as
    # historical data and are NOT re-run.
    "V3Fast":    {"env": "PIVOTER_RUN_REGION_V3FAST"},
    # V3Fast_NP: same optimized path as V3Fast but with Private Cloud OFF —
    # strict ablation study for the Private Cloud optimization. Renamed
    # from "V3_NP" so the new optimized ablation data is separate from the
    # legacy V3_NP data that was collected on the pre-optimization code.
    "V3Fast_NP": {"env": "PIVOTER_RUN_REGION_V3FAST", "extra": {"PIVOTER_V3_NO_PRIVATE": "1"}},
    # V3Fast_NoCPI: ablation for the CPI support formula (Theorem 6.1).
    # Identical to V3Fast except Step 4 (initial support) enumerates every
    # s-subset of every region instead of using the CPI convolution.
    # Expected to time out on graphs with large MC size × large s; that's
    # the intended evidence of the CPI formula's value.
    "V3Fast_NoCPI": {"env": "PIVOTER_RUN_REGION_V3NOCPI"},
}

# ============ Helpers ============
def link_graphs():
    os.makedirs("graphs", exist_ok=True)
    missing = []
    for g in GRAPHS:
        f = f"graphs/{g}.edges"
        src = f"{DATADIR}/{g}.edges"
        if os.path.exists(f):
            continue
        if os.path.exists(src):
            os.symlink(src, f)
            print(f"  Linked {g}.edges")
        else:
            missing.append((g, f, src))
    if missing:
        print("\n!! Missing graph files — these entries will NOT run:")
        for g, f, src in missing:
            print(f"    {g}: neither {f} nor {src} exists")
        print("  Fix: copy/download the .edges file, or remove the entry from SERVER_GRAPHS.\n")

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

def get_proc_rss_gb(pid):
    """Get RSS of a process in GB from /proc/[pid]/status."""
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1]) / 1024 / 1024  # kB → GB
    except:
        pass
    return 0

def get_used_mem_gb():
    """Get used memory in GB from /proc/meminfo."""
    try:
        with open("/proc/meminfo") as f:
            info = {}
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    info[parts[0].rstrip(":")] = int(parts[1])
            total = info.get("MemTotal", 0)
            avail = info.get("MemAvailable", 0)
            return (total - avail) / 1024 / 1024
    except:
        return 0

def probe_max_clique(graph):
    """Compute max clique size with NetworkX (authoritative, no C++ probe).
    Returns 0 on failure so caller can warn rather than cache a wrong value."""
    gf = f"graphs/{graph}.edges"
    if not os.path.exists(gf):
        return 0
    try:
        import networkx as nx
    except ImportError:
        print(f"  !! networkx not installed — run: pip install networkx")
        return 0
    try:
        G = nx.Graph()
        with open(gf) as f:
            header = f.readline().split()
            # skip the "n m" header; add edges from the rest
            for ln in f:
                parts = ln.split()
                if len(parts) < 2:
                    continue
                u, v = int(parts[0]), int(parts[1])
                if u != v:
                    G.add_edge(u, v)
        best = 0
        for c in nx.find_cliques(G):
            if len(c) > best:
                best = len(c)
        return best
    except Exception as e:
        print(f"  !! networkx probe failed for {graph}: {e}")
        return 0

def load_existing():
    """Load already-completed results to skip."""
    done = set()
    if os.path.exists(OUTCSV):
        with open(OUTCSV) as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    key = (row["graph"], int(row["r"]), int(row["s"]), row["algo"])
                    done.add(key)
                except (ValueError, KeyError):
                    pass
    return done

def extract_timing(txt):
    """Extract total_ms, peel_ms, mem_kB from log text."""
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', txt)
    m_peel = re.search(r'Peeling time:\s*([\d.]+)', txt)
    m_mem = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)\s*kB', txt)
    total_ms = float(m_total.group(1)) if m_total else -1.0
    peel_ms = float(m_peel.group(1)) if m_peel else -1.0
    mem_kB = float(m_mem.group(1)) if m_mem else -1.0
    return total_ms, peel_ms, mem_kB

def write_result(graph, r, s, algo, status, wall_ms=-1, total_ms=-1, peel_ms=-1, mem_kB=-1):
    """Append one result row to CSV."""
    with open(OUTCSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writerow({
            "graph": graph, "r": r, "s": s, "algo": algo, "status": status,
            "wall_ms": f"{wall_ms:.1f}" if wall_ms >= 0 else "",
            "total_ms": f"{total_ms:.1f}" if total_ms >= 0 else "",
            "peel_ms": f"{peel_ms:.1f}" if peel_ms >= 0 else "",
            "mem_kB": f"{mem_kB:.0f}" if mem_kB >= 0 else "",
        })

FIELDNAMES = ["graph", "r", "s", "algo", "status", "wall_ms", "total_ms", "peel_ms", "mem_kB"]

# ============ Main ============
def main():
    print("=" * 60)
    print("  V3 Full Benchmark")
    print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    link_graphs()
    build()
    LOGDIR.mkdir(exist_ok=True)

    # Setup CSV header
    if not os.path.exists(OUTCSV):
        with open(OUTCSV, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=FIELDNAMES).writeheader()

    # Probe max clique sizes (cached to JSON)
    cache_file = Path("bench_v3_max_cliques.json")
    max_cliques = {}
    if cache_file.exists():
        max_cliques = json.loads(cache_file.read_text())
        print(f"\nLoaded max clique cache: {cache_file}")

        # Invalidate cache entries whose .edges file is newer than the cache
        # (protects against stale values after edits like dedup / regeneration).
        cache_mtime = cache_file.stat().st_mtime
        stale = []
        for g in list(max_cliques.keys()):
            p = f"graphs/{g}.edges"
            if os.path.exists(p) and os.path.getmtime(p) > cache_mtime:
                stale.append(g)
                del max_cliques[g]
        if stale:
            print(f"  Invalidated stale cache entries (file newer than cache): {stale}")

        # 需要 probe 的图 = GRAPHS 里当前没有 cache 值的
        missing_graphs = [g for g in GRAPHS if g not in max_cliques and os.path.exists(f"graphs/{g}.edges")]
        if missing_graphs:
            print(f"\nProbing missing graphs (incremental): {missing_graphs}")
            with ProcessPoolExecutor(max_workers=len(missing_graphs)) as ex:
                futures = {ex.submit(probe_max_clique, g): g for g in missing_graphs}
                for f in futures:
                    g = futures[f]
                    mc = f.result()
                    if mc == 0:
                        print(f"  !! probe failed for {g}: no maxSize in output — NOT caching.")
                        print(f"     Run manually to see the real error:")
                        print(f"       PIVOTER_RUN_REGION_V3=1 PIVOTER_V3_NO_PRIVATE=1 {BIN} graphs/{g}.edges 3 4 degen")
                    else:
                        max_cliques[g] = mc
            # 更新并重新保存缓存
            cache_file.write_text(json.dumps(max_cliques, indent=2))
            print(f"  Updated cache saved to {cache_file}")

        for g in GRAPHS: # 按 GRAPHS 顺序打印，明确哪些被跳过
            if g in max_cliques:
                mc = max_cliques[g]
                n_combos = sum(s - 3 for s in range(4, mc + 1))
                print(f"  {g}: max_clique={mc}, jobs={n_combos * len(ALGOS)}")
            elif not os.path.exists(f"graphs/{g}.edges"):
                print(f"  {g}: SKIPPED (file missing)")
            else:
                print(f"  {g}: SKIPPED (not in cache, probe may have failed)")
    else:
        print("\nProbing max clique sizes (parallel)...")
        graphs_to_probe = [g for g in GRAPHS if os.path.exists(f"graphs/{g}.edges")]
        if not graphs_to_probe:
            print("  No graphs found! Check DATADIR and graph files.")
            sys.exit(1)
        with ProcessPoolExecutor(max_workers=len(graphs_to_probe)) as ex:
            futures = {ex.submit(probe_max_clique, g): g for g in graphs_to_probe}
            for f in futures:
                g = futures[f]
                mc = f.result()
                if mc == 0:
                    print(f"  !! probe failed for {g}: no maxSize in output — NOT caching.")
                    continue
                max_cliques[g] = mc
                n_combos = sum(s - 3 for s in range(4, mc + 1))
                print(f"  {g}: max_clique={mc}, jobs={n_combos * len(ALGOS)}")
        cache_file.write_text(json.dumps(max_cliques, indent=2))
        print(f"  Saved to {cache_file}")

    done = load_existing()
    total_jobs = sum(
        sum(s - 3 for s in range(4, max_cliques[g] + 1)) * len(ALGOS)
        for g in max_cliques
    )
    print(f"\nTotal jobs: {total_jobs}, already done: {len(done)}")
    print(f"Max workers: {MAX_WORKERS}, mem limit: {MEM_LIMIT_GB}GB, "
          f"kill: {MEM_KILL_GB}GB, timeout: {TIMEOUT}s")

    # Timeout tracking: (graph, algo, s) → min r that timed out
    timeout_at = defaultdict(lambda: float('inf'))

    # Load timeout info from existing results and re-apply propagation
    if os.path.exists(OUTCSV):
        # First pass: collect direct timeouts
        direct_timeouts = []  # (graph, algo, s, r)
        with open(OUTCSV) as f:
            for row in csv.DictReader(f):
                if row.get("status") in ("TIMEOUT", "OOM"):
                    try:
                        direct_timeouts.append((
                            row["graph"], row["algo"],
                            int(row["s"]), int(row["r"])
                        ))
                    except:
                        pass
        # Apply with propagation (ST/REF propagate across all s for same r)
        for g, an, ss, rr in direct_timeouts:
            if an in ("ST", "REF"):
                for sf in range(4, max_cliques.get(g, 0) + 1):
                    timeout_at[(g, an, sf)] = min(timeout_at[(g, an, sf)], rr)
            else:
                timeout_at[(g, an, ss)] = min(timeout_at[(g, an, ss)], rr)
        n_skip = sum(1 for k, v in timeout_at.items() if v < float('inf'))
        print(f"  Loaded {len(direct_timeouts)} timeouts, {n_skip} skip rules")

    # Track running processes: [(Popen, graph, r, s, algo, start_time)]
    running = []
    launched = 0
    skipped = 0
    retry_count = defaultdict(int)  # (graph, r, s, algo) → how many times re-queued
    MAX_RETRIES = 2

    # Graceful shutdown
    shutdown = False
    def handle_signal(sig, frame):
        nonlocal shutdown
        print("\n[SIGINT] Shutting down — killing all children...")
        shutdown = True
        for proc, *_ in running:
            try: proc.kill()
            except: pass
    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)

    def reap():
        """Check for finished processes, record results."""
        nonlocal launched
        still = []
        for proc, g, rr, ss, an, t0 in running:
            ret = proc.poll()
            if ret is None:
                still.append((proc, g, rr, ss, an, t0))
                continue
            try:
                txt, _ = proc.communicate(timeout=5)
            except:
                txt = ""
            if ret == -9 or ret == 137:
                status = "OOM"
            elif ret == 0:
                status = "OK"
            else:
                # SIGSEGV=-11, SIGBUS=-10, SIGFPE=-8 etc
                sig_names = {-11: "SIGSEGV", -10: "SIGBUS", -8: "SIGFPE", -6: "SIGABRT", -4: "SIGILL"}
                status = f"ERROR({ret}{'/' + sig_names[ret] if ret in sig_names else ''})"
            (LOGDIR / f"{g}_r{rr}_s{ss}_{an}.log").write_text(txt)
            wall_ms = (time.time() - t0) * 1000
            total_ms, peel_ms, mem_kB = extract_timing(txt)
            t_str = f"{wall_ms:.0f}ms" if wall_ms >= 0 else "N/A"
            m_str = f"{mem_kB/1024:.0f}MB" if mem_kB >= 0 else ""
            print(f"  {an:>6} {g} r={rr} s={ss} {status} wall={t_str} {m_str}", flush=True)
            if status.startswith("ERROR") and txt:
                # Print last few lines of stdout/stderr to help diagnose
                tail = "\n".join(txt.strip().splitlines()[-5:])
                if tail:
                    print(f"         last output: {tail}", flush=True)
            write_result(g, rr, ss, an, status, wall_ms, total_ms, peel_ms, mem_kB)
            done.add((g, rr, ss, an))
            launched += 1
            if status in ("TIMEOUT", "OOM"):
                propagate_timeout(g, an, ss, rr)
        running.clear()
        running.extend(still)

    def check_timeouts():
        """Kill processes exceeding TIMEOUT."""
        now = time.time()
        for i in range(len(running) - 1, -1, -1):
            proc, g, rr, ss, an, t0 = running[i]
            if now - t0 > TIMEOUT + 10:
                try:
                    proc.kill()
                    txt, _ = proc.communicate(timeout=10)
                except:
                    txt = ""
                    try: proc.wait(timeout=5)
                    except: pass
                print(f"  {an:>6} {g} r={rr} s={ss} TIMEOUT", flush=True)
                (LOGDIR / f"{g}_r{rr}_s{ss}_{an}.log").write_text(
                    txt + "\nTIMEOUT by scheduler")
                write_result(g, rr, ss, an, "TIMEOUT")
                done.add((g, rr, ss, an))
                propagate_timeout(g, an, ss, rr)
                running.pop(i)

    def check_proc_mem():
        """Kill any process whose RSS exceeds PER_PROC_MEM_GB."""
        for i in range(len(running) - 1, -1, -1):
            proc, g, rr, ss, an, t0 = running[i]
            rss = get_proc_rss_gb(proc.pid)
            if rss > PER_PROC_MEM_GB:
                try:
                    proc.kill()
                    proc.communicate(timeout=10)
                except:
                    try: proc.wait(timeout=5)
                    except: pass
                print(f"  [OOM] {an} {g} r={rr} s={ss} RSS={rss:.0f}GB > {PER_PROC_MEM_GB}GB", flush=True)
                (LOGDIR / f"{g}_r{rr}_s{ss}_{an}.log").write_text(f"OOM: RSS={rss:.0f}GB")
                write_result(g, rr, ss, an, "OOM")
                done.add((g, rr, ss, an))
                propagate_timeout(g, an, ss, rr)
                running.pop(i)

    def kill_newest():
        """Kill most recent process for OOM protection."""
        if not running:
            return
        proc, g, rr, ss, an, t0 = running.pop()
        try:
            proc.kill()
            proc.communicate(timeout=10)
        except:
            try: proc.wait(timeout=5)
            except: pass
        mem = get_used_mem_gb()
        print(f"  [KILL] {an} {g} r={rr} s={ss} (mem={mem:.0f}GB)", flush=True)
        if len(running) == 0:
            # Only job → genuine OOM
            write_result(g, rr, ss, an, "OOM")
            done.add((g, rr, ss, an))
            propagate_timeout(g, an, ss, rr)
        else:
            # Multiple jobs → re-queue (with retry limit)
            key = (g, rr, ss, an)
            retry_count[key] += 1
            if retry_count[key] > MAX_RETRIES:
                write_result(g, rr, ss, an, "OOM")
                done.add(key)
                propagate_timeout(g, an, ss, rr)
                print(f"    → OOM (retried {MAX_RETRIES} times)", flush=True)
            else:
                retry_queue.append(key)
                print(f"    → re-queued ({retry_count[key]}/{MAX_RETRIES})", flush=True)

    def should_skip(g, rr, ss, an):
        return rr >= timeout_at[(g, an, ss)]

    def launch(g, rr, ss, an):
        env = {**os.environ, ALGOS[an]["env"]: "1"}
        if "extra" in ALGOS[an]:
            env.update(ALGOS[an]["extra"])
        proc = subprocess.Popen(
            [BIN, f"graphs/{g}.edges", str(rr), str(ss)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env
        )
        running.append((proc, g, rr, ss, an, time.time()))

    def propagate_timeout(g, an, ss, rr):
        """ST: propagate along s for same r (cliqueIndex depends on r, not s).
        V3/V3_NP: no propagation (performance depends on both r and s)."""
        if an in ("ST", "REF"):
            # cliqueIndex is O(C(n,r)) — if r=K fails, it fails for ALL s
            for sf in range(4, max_cliques.get(g, 0) + 1):
                timeout_at[(g, an, sf)] = min(timeout_at[(g, an, sf)], rr)

    # ---- Build global job queue: ordered by (graph, s, r, algo) ----
    # This order ensures timeout propagation works (smaller s before larger s)
    retry_queue = []
    all_jobs = []
    for graph in GRAPHS:
        if graph not in max_cliques:
            continue
        max_k = max_cliques[graph]
        for s in range(4, max_k + 1):
            for r in range(3, s):
                for algo_name in ALGOS:
                    all_jobs.append((graph, r, s, algo_name))

    ji = 0
    while (ji < len(all_jobs) or retry_queue or running) and not shutdown:
        reap()
        check_timeouts(); check_proc_mem()

        # OOM protection
        while get_used_mem_gb() > MEM_KILL_GB and running:
            kill_newest()
            time.sleep(2)

        # Pick next job: retry queue first, then main queue
        job = None
        if retry_queue:
            job = retry_queue.pop(0)
        elif ji < len(all_jobs):
            job = all_jobs[ji]
            ji += 1
        else:
            # No more jobs to launch, just wait for running ones
            if running:
                time.sleep(POLL_SEC)
            continue

        g, rr, ss, an = job

        # Skip if already done
        if (g, rr, ss, an) in done:
            continue

        # Skip if timed out
        if should_skip(g, rr, ss, an):
            skipped += 1
            write_result(g, rr, ss, an, "SKIP_TIMEOUT")
            done.add((g, rr, ss, an))
            continue

        # Wait for capacity: workers < MAX and memory < limit
        while (len(running) >= MAX_WORKERS or get_used_mem_gb() >= MEM_LIMIT_GB) and not shutdown:
            reap()
            check_timeouts(); check_proc_mem()
            while get_used_mem_gb() > MEM_KILL_GB and running:
                kill_newest()
                time.sleep(2)
            if len(running) >= MAX_WORKERS or get_used_mem_gb() >= MEM_LIMIT_GB:
                time.sleep(POLL_SEC)

        if shutdown:
            break

        # Re-check skip after waiting (timeout_at may have been updated)
        if (g, rr, ss, an) in done or should_skip(g, rr, ss, an):
            if (g, rr, ss, an) not in done:
                skipped += 1
                write_result(g, rr, ss, an, "SKIP_TIMEOUT")
                done.add((g, rr, ss, an))
            continue

        launch(g, rr, ss, an)
        # Dynamic throttle: brief pause, longer if memory is climbing
        mem_now = get_used_mem_gb()
        if mem_now > MEM_LIMIT_GB * 0.8:    # >240GB: slow down
            time.sleep(3)
        elif mem_now > MEM_LIMIT_GB * 0.5:   # >150GB: moderate
            time.sleep(1)
        else:
            time.sleep(SETTLE_SEC)            # low mem: fast launch

    # Drain remaining
    while running and not shutdown:
        reap()
        check_timeouts(); check_proc_mem()
        time.sleep(POLL_SEC)
    reap()

    print(f"\n{'='*60}")
    print(f"  DONE. Launched: {launched}, skipped: {skipped}")
    print(f"  Results: {OUTCSV}")
    print(f"  Logs: {LOGDIR}/")
    print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()