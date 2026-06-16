#!/usr/bin/env python3
"""
bench_arb.py — ARB-NUCLEUS-DECOMP parallel baseline sweep (single-thread,
both variants) for the RegND* paper comparison.

Runs the GBBS NucleusDecomposition_main binaries on the SIX current paper
graphs, converted to GBBS adjacency format by edges_to_gbbs_adj.py (we do
NOT reuse the stale /data/wenqianz/*.adj which are a different dataset).

Two variants, matching our two algorithms:
  ARB     (~/arb-nucleus-decomp)    core numbers only      <-> RegND*
  ARB-Hi  (~/arb-nucleus-hierarchy) with hierarchy output  <-> RegND-H

Each cell runs SINGLE-THREADED (PARLAY_NUM_THREADS=1) so the numbers align
with the paper's single-thread RegND*/CND numbers. The sweep runs several
cells concurrently to finish faster; --workers controls concurrency and a
memory gate avoids over-subscription. NOTE: co-running cells share memory
bandwidth, which can inflate timings by a few percent; for the final
figure cells, a sequential re-measure (--workers 1 on a cell whitelist)
is the clean-number path, mirroring the hierarchy A/B re-measure.

Each run is wrapped with /usr/bin/time -v; we record the algorithm's
internal "### Running Time", the wrapper's elapsed wall time, user time,
peak RSS, and exit/signal. Resumable: cells already terminal in the CSV
are skipped. Skip-floor: once (graph, variant, r) times out or OOMs at
some s, higher s for that (graph, variant, r) are skipped.

Usage:
  nohup python3 bench_arb.py > bench_arb.log 2>&1 &
  BENCH_ARB_WORKERS=12 BENCH_ARB_RGRID=3,4,5,6,7 python3 bench_arb.py
"""
from __future__ import annotations
import csv, os, re, signal, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HOME = Path.home()
BINS = {
    "ARB":    HOME / "arb-nucleus-decomp/bazel-bin/benchmarks/NucleusDecomposition/NucleusDecomposition_main",
    "ARB-Hi": HOME / "arb-nucleus-hierarchy/bazel-bin/benchmarks/NucleusDecomposition/NucleusDecomposition_main",
}
ADJ_DIR  = Path("/data/wenqianz/arb_adj")          # our-graph .adj files live here
GRAPHS   = ["dblp-core30", "ca-GrQc", "ca-HepPh", "ca-CondMat", "com-dblp", "web-it-2004"]
RGRID    = [int(x) for x in os.getenv("BENCH_ARB_RGRID", "3,4,5,6,7").split(",")]
S_CAP    = int(os.getenv("BENCH_ARB_SCAP", "50"))  # safety ceiling on s
TIMEOUT  = int(os.getenv("BENCH_ARB_TIMEOUT", "3600"))
WORKERS  = int(os.getenv("BENCH_ARB_WORKERS", "24"))
MIN_FREE_GB = int(os.getenv("BENCH_ARB_MIN_FREE_GB", "80"))  # memory gate
# RLIMIT_AS per cell (GB). ARB-Hi's hierarchy/connectivity structure can
# balloon to 100-200GB even on tiny dense graphs, so cap each cell at the
# paper's 250GB per-run budget; a cell exceeding it is recorded OOM (a
# valid "out of budget" result) instead of taking down the 503GB box.
MEMCAP_GB = int(os.getenv("BENCH_ARB_MEMCAP_GB", "0"))  # 0 = no cap
# restrict which variants to run (e.g. ARB-Hi alone after ARB is done)
VWANT = set(v.strip() for v in os.getenv("BENCH_ARB_VARIANTS", "ARB,ARB-Hi").split(","))
OUTCSV   = Path(os.getenv("BENCH_ARB_OUTCSV", "paper_data/bench_arb.csv"))
LOGDIR   = Path("bench_arb_logs")
FIELDS   = ["graph", "r", "s", "variant", "status",
            "alg_s", "user_s", "wall_s", "max_rss_mb", "exit_status"]

ALG_RE  = re.compile(r"### Running Time:\s*([0-9.eE+-]+)")
USER_RE = re.compile(r"User time \(seconds\):\s*([0-9.]+)")
WALL_RE = re.compile(r"Elapsed \(wall clock\) time .*?:\s*([0-9:.]+)")
RSS_RE  = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
EXIT_RE = re.compile(r"Exit status:\s*(\d+)")
SIG_RE  = re.compile(r"Command terminated by signal\s+(\d+)")


def wall_to_sec(s: str) -> float:
    parts = s.split(":")
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        return float(parts[0])
    except ValueError:
        return -1.0


def load_done(path: Path) -> set:
    done = set()
    if path.exists():
        with path.open() as f:
            for row in csv.DictReader(f):
                done.add((row["graph"], int(row["r"]), int(row["s"]), row["variant"]))
    return done


def omega_of(graph: str) -> int:
    """Max s to attempt = clique number, read from the converted .adj is
    not possible; use a generous per-graph ceiling and rely on timeout-
    skip. S_CAP bounds it globally."""
    return S_CAP


def free_gb() -> float:
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024.0 * 1024.0)
    except Exception:
        pass
    return 1e9   # no /proc/meminfo: don't gate


def mem_gate():
    """Block until at least MIN_FREE_GB is available; throttles heavy
    cells (web-it, com-dblp at high r) without starving the many light
    cells on the small graphs."""
    waited = 0
    while free_gb() < MIN_FREE_GB:
        time.sleep(5); waited += 5
        if waited and waited % 300 == 0:
            print(f"  [mem-gate] waiting, free={free_gb():.0f}GB "
                  f"< {MIN_FREE_GB}GB ({waited}s)", flush=True)


def run_cell(graph: str, variant: str, r: int, s: int):
    binpath = BINS[variant]
    adj = ADJ_DIR / f"{graph}.adj"
    LOGDIR.mkdir(parents=True, exist_ok=True)
    mem_gate()
    log = LOGDIR / f"{graph}_{variant}_r{r}_s{s}.timev.log"
    # README canonical invocation: -r/-ss select the clique sizes, -relabel
    # + -efficient 1 + -tt 5 select the actual decomposition path. WITHOUT
    # -tt (default tt=0) the binary parses args and exits in ~ms WITHOUT
    # decomposing -- the first sweep was entirely no-ops (caught because
    # 4ms on a graph where CND/RegND* take minutes is implausible; the
    # correct-flag max core 236 on ca-HepPh matches RegND*).
    cmd = ["/usr/bin/time", "-v", str(binpath),
           "-rounds", "1", "-s", "-r", str(r), "-ss", str(s),
           "-relabel", "-efficient", "1", "-tt", "5", str(adj)]
    env = dict(os.environ, PARLAY_NUM_THREADS="1")

    def child_setup():
        os.setsid()
        if MEMCAP_GB > 0:
            import resource
            cap = MEMCAP_GB * 1024 * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    t0 = time.time()
    try:
        # GBBS prints "### Running Time" to STDOUT; /usr/bin/time -v prints
        # User/wall/maxRSS to STDERR. Capture both.
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, env=env,
                             preexec_fn=child_setup, text=True)
        try:
            out, err = p.communicate(timeout=TIMEOUT)
            status = "OK"
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgid(p.pid), signal.SIGKILL)
            out, err = p.communicate()
            status = "TIMEOUT"
    except Exception as e:
        return dict(graph=graph, r=r, s=s, variant=variant, status="ERROR",
                    alg_s=-1, user_s=-1, wall_s=time.time() - t0,
                    max_rss_mb=-1, exit_status=str(e)[:40])
    log.write_text(out + "\n----TIME-V----\n" + err)
    alg  = ALG_RE.search(out)
    user = USER_RE.search(err)
    wall = WALL_RE.search(err)
    rss  = RSS_RE.search(err)
    ex   = EXIT_RE.search(err)
    sig  = SIG_RE.search(err)
    exit_status = int(ex.group(1)) if ex else 0
    if sig:
        # SIGKILL from our timeout already labeled; other signals = OOM/crash.
        # RLIMIT_AS exhaustion surfaces as std::bad_alloc -> SIGABRT(6), or a
        # null-deref SIGSEGV(11), or kernel SIGKILL(9).
        if status != "TIMEOUT":
            status = "OOM" if sig.group(1) in ("6", "9", "11") else "ERROR"
        exit_status = 128 + int(sig.group(1))
    elif status == "OK" and ("bad_alloc" in out or "bad_alloc" in err
                             or "out of memory" in (out + err).lower()):
        status = "OOM"
    elif status == "OK" and not alg:
        status = "ERROR"   # finished but no timing line: treat as failure
    return dict(graph=graph, r=r, s=s, variant=variant, status=status,
                alg_s=float(alg.group(1)) if alg else -1,
                user_s=float(user.group(1)) if user else -1,
                wall_s=wall_to_sec(wall.group(1)) if wall else (time.time() - t0),
                max_rss_mb=(int(rss.group(1)) / 1024.0) if rss else -1,
                exit_status=exit_status)


def main():
    for v, b in BINS.items():
        if not b.exists():
            sys.exit(f"missing {v} binary: {b}")
    OUTCSV.parent.mkdir(parents=True, exist_ok=True)
    done = load_done(OUTCSV)
    new_file = not OUTCSV.exists()
    fout = OUTCSV.open("a", newline="")
    w = csv.DictWriter(fout, fieldnames=FIELDS)
    if new_file:
        w.writeheader(); fout.flush()

    # build job list per (graph, variant, r) chain so skip-floor is local
    chains = []
    for graph in GRAPHS:
        if not (ADJ_DIR / f"{graph}.adj").exists():
            print(f"[skip] no .adj for {graph}", flush=True); continue
        for variant in BINS:
            if variant not in VWANT:
                continue
            for r in RGRID:
                chains.append((graph, variant, r))

    print(f"[arb] {len(chains)} (graph,variant,r) chains, s up to {S_CAP}, "
          f"{WORKERS} workers, timeout {TIMEOUT}s, single-thread", flush=True)

    import threading
    wlock = threading.Lock()
    counter = {"n": 0}

    def emit(res):
        # Crash-safe: write each cell immediately under a lock so a long
        # slow chain (ARB-Hi) is never lost on kill, and progress is
        # visible live. Per-chain batching delayed all rows until the
        # whole (slow) chain returned.
        with wlock:
            w.writerow(res); fout.flush(); counter["n"] += 1
            print(f"  {res['graph']:14} {res['variant']:7} "
                  f"r={res['r']} s={res['s']:2} {res['status']:8} "
                  f"alg={res['alg_s']}s wall={res['wall_s']:.1f}s "
                  f"rss={res['max_rss_mb']:.0f}MB", flush=True)

    # each chain is sequential internally (skip-floor); chains run concurrently
    def run_chain(graph, variant, r):
        for s in range(r + 1, omega_of(graph) + 1):
            if (graph, r, s, variant) in done:
                continue
            res = run_cell(graph, variant, r, s)
            emit(res)
            if res["status"] in ("TIMEOUT", "OOM"):
                break   # skip-floor: stop this chain at the first death

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(run_chain, g, v, r) for g, v, r in chains]
        for fut in as_completed(futs):
            fut.result()   # surface exceptions; rows already emitted
    fout.close()
    print(f"[arb] {counter['n']} cells written to {OUTCSV} in "
          f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
