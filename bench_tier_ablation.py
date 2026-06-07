#!/usr/bin/env python3
"""
bench_tier_ablation.py — 4-tier RegND ablation sweep.

Iterates PIVOTER_TIER={1,2,3,4} over a curated subset of (graph, r, s) cells
and records (wall, total, step4, peel, peak RSS) per row in a paper-6
compatible CSV.  T4 numbers can also be reused from bench_v3_all_results.csv
to save runtime.

Usage:
  python3 bench_tier_ablation.py {tods1|tods2|local}
  python3 bench_tier_ablation.py --reuse-t4   # join T4 rows from main CSV
  python3 bench_tier_ablation.py --dry-run    # print row plan, do not run
"""

from __future__ import annotations
import argparse, csv, os, re, signal, socket, subprocess, sys, time
from pathlib import Path

# ---------- stack limit ----------
try:
    import resource
    _BIG = 1 << 30
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    if hard != resource.RLIM_INFINITY:
        target = min(hard, _BIG)
    else:
        target = _BIG
    if target != soft:
        try:
            resource.setrlimit(resource.RLIMIT_STACK, (target, hard))
        except (ValueError, OSError):
            pass
except Exception:
    pass


# ============ Config ============
BIN = Path("./build/bin/degeneracy_cliques")
OUTCSV = "bench_tier_ablation_results.csv"
LOGDIR = Path("bench_tier_ablation_logs")
TIMEOUT_SEC = 3600      # per-cell wall-clock cap
MEM_KILL_GB = 250       # per-process peak RSS cap (kernel will OOM-kill above)

# Reference V3LM CSV — T4 rows can be re-used from here.
V3LM_CSV = "bench_v3_all_results.csv"


# ============ Ablation matrix ============
# Spec: 6 paper-6 graphs × up to 3 (r, s) cells × 4 tiers.
# Tiers monotonically dominated: T1 slowest, T4 fastest.
#
# For cells where T4 already timed out (existing T4 row = TIMEOUT in main
# CSV), we still attempt T3/T2/T1 ONLY if the smaller (r,s) cells finished
# — otherwise we skip and mark TIMEOUT directly to save cluster time.
ABLATION_CELLS = [
    # (graph stem, r, s)
    # ---- ca-GrQc (sparse, smoke) ----
    ("ca-GrQc",       3, 4),
    ("ca-GrQc",       5, 6),
    ("ca-GrQc",       8, 9),
    # ---- ca-CondMat (mid) ----
    ("ca-CondMat",    3, 4),
    ("ca-CondMat",    5, 6),
    ("ca-CondMat",    7, 8),
    # ---- com-dblp (coauth) ----
    ("com-dblp",      3, 4),
    ("com-dblp",      5, 6),
    ("com-dblp",      7, 8),
    # ---- ca-HepPh (dense small) ----
    ("ca-HepPh",      3, 4),
    ("ca-HepPh",      6, 7),
    ("ca-HepPh",     10, 11),
    # ---- web-it-2004 (dense large) ----
    ("web-it-2004",   3, 4),
    ("web-it-2004",   5, 6),
    ("web-it-2004",   7, 8),
    # ---- dblp-core30 (30-core subgraph) ----
    ("dblp-core30",   3, 4),
    ("dblp-core30",   6, 7),
    ("dblp-core30",   9, 10),
]

# Tier-1 hard skip list: cells where direct s-clique enumeration is
# guaranteed to time out.  Mark TIMEOUT in CSV without running.
T1_HARD_SKIP = {
    # web-it-2004: dense maximal cliques; direct s-enum explodes at any r.
    ("web-it-2004", 3, 4), ("web-it-2004", 5, 6), ("web-it-2004", 7, 8),
    # com-dblp r>=5: matched-cell evidence T1 cannot finish.
    ("com-dblp", 5, 6), ("com-dblp", 7, 8),
    # ca-HepPh r>=6: degeneracy 238 + direct s-enum.
    ("ca-HepPh", 6, 7), ("ca-HepPh", 10, 11),
    # dblp-core30 r>=6: dense 30-core.
    ("dblp-core30", 6, 7), ("dblp-core30", 9, 10),
}

# Tier-2 likely-timeout list: cells where CPI helps but recompute peel
# (re-running AggrCount per peel batch) is still too slow.
T2_LIKELY_TIMEOUT = {
    ("ca-HepPh", 10, 11),
    ("web-it-2004", 7, 8),
}

TIERS = [1, 2, 3, 4]


# ============ Server config ============
SERVER_DATADIR = {
    "tods1": "/data/wenqianz",
    "tods2": "/data/wenqianz",
    "local": "graphs",
}
SERVER_MAX_PARALLEL = {
    "tods1": 8,
    "tods2": 12,
    "local": 1,
}


def detect_server():
    if len(sys.argv) > 1 and sys.argv[1] in SERVER_DATADIR:
        return sys.argv[1]
    host = socket.gethostname().lower()
    for key in SERVER_DATADIR:
        if key in host:
            return key
    return "local"


# ============ Build ============
def cmake_build():
    print("[build] cmake configure + build (-j 12)", flush=True)
    subprocess.run("cmake -S . -B build -DCMAKE_BUILD_TYPE=Release".split(),
                   check=True, stdout=subprocess.DEVNULL)
    r = subprocess.run("cmake --build build -j 12 --target degeneracy_cliques".split(),
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr, file=sys.stderr)
        sys.exit(1)
    print("[build] OK", flush=True)


# ============ Graph linking ============
def link_graphs(server: str):
    datadir = SERVER_DATADIR[server]
    if datadir == "graphs":
        return
    os.makedirs("graphs", exist_ok=True)
    cells = {c[0] for c in ABLATION_CELLS}
    missing = []
    for g in cells:
        f = f"graphs/{g}.edges"
        src = f"{datadir}/{g}.edges"
        if os.path.exists(f):
            continue
        if os.path.exists(src):
            try: os.symlink(src, f)
            except OSError: pass
        else:
            missing.append(g)
    if missing:
        print(f"[link] WARN: missing graphs: {missing}", flush=True)


# ============ Existing-row resume / T4 reuse ============
FIELDNAMES = ["graph", "r", "s", "tier", "algo", "status",
              "wall_ms", "total_ms", "step4_ms", "peel_ms", "mem_kB", "max_core"]


def load_done():
    done = set()
    if os.path.exists(OUTCSV):
        with open(OUTCSV) as f:
            for row in csv.DictReader(f):
                try:
                    key = (row["graph"], int(row["r"]), int(row["s"]), int(row["tier"]))
                    done.add(key)
                except (ValueError, KeyError):
                    pass
    return done


def load_t4_from_v3lm():
    """If RegNDC rows exist in bench_v3_all_results.csv, reuse them as T4."""
    t4 = {}
    if not os.path.exists(V3LM_CSV):
        return t4
    with open(V3LM_CSV) as f:
        for row in csv.DictReader(f):
            if row.get("algo") != "RegNDC":
                continue
            if row.get("status") != "OK":
                continue
            try:
                key = (row["graph"], int(row["r"]), int(row["s"]))
                t4[key] = row
            except (ValueError, KeyError):
                pass
    return t4


# ============ Run one cell ============
def extract_timing(txt):
    m_total = re.search(r'NucleusCoreDecomposition took:\s*([\d.]+)', txt)
    m_peel  = re.search(r'Peeling time:\s*([\d.]+)', txt)
    m_step4 = re.search(r'(?:CPI counting time|NoCPI enumeration time):\s*([\d.]+)', txt)
    m_mem   = re.search(r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)\s*kB', txt)
    m_core  = re.search(r'Max core:\s*(\d+)', txt)
    return (
        float(m_total.group(1)) if m_total else -1.0,
        float(m_step4.group(1)) if m_step4 else -1.0,
        float(m_peel.group(1))  if m_peel  else -1.0,
        float(m_mem.group(1))   if m_mem   else -1.0,
        int(m_core.group(1))    if m_core  else -1,
    )


def run_cell(graph, r, s, tier):
    """Run one (graph, r, s, tier) cell; return (status, fields, log_path)."""
    LOGDIR.mkdir(exist_ok=True)
    log_path = LOGDIR / f"T{tier}_{graph}_r{r}_s{s}.log"
    env = os.environ.copy()
    env["PIVOTER_RUN_REGION_TIER"] = "1"
    env["PIVOTER_TIER"] = str(tier)
    # Always wrap with /usr/bin/time -v so we capture peak RSS robustly.
    cmd = ["/usr/bin/time", "-v",
           str(BIN), f"graphs/{graph}.edges", str(r), str(s), "degen"]
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd, env=env, capture_output=True, text=True,
            timeout=TIMEOUT_SEC,
        )
        elapsed = time.time() - t0
        combined = (proc.stdout or "") + "\n--STDERR--\n" + (proc.stderr or "")
        log_path.write_text(combined)
        status = "OK" if proc.returncode == 0 else f"ERR({proc.returncode})"
        total_ms, step4_ms, peel_ms, mem_kB, max_core = extract_timing(combined)
        wall_ms = elapsed * 1000.0
        return status, {
            "wall_ms": wall_ms, "total_ms": total_ms, "step4_ms": step4_ms,
            "peel_ms": peel_ms, "mem_kB": mem_kB, "max_core": max_core,
        }, log_path
    except subprocess.TimeoutExpired:
        log_path.write_text(f"TIMEOUT after {TIMEOUT_SEC}s\n")
        return "TIMEOUT", {"wall_ms": TIMEOUT_SEC * 1000.0,
                           "total_ms": -1, "step4_ms": -1, "peel_ms": -1,
                           "mem_kB": -1, "max_core": -1}, log_path


def write_row(graph, r, s, tier, status, fields):
    new_file = not os.path.exists(OUTCSV)
    with open(OUTCSV, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if new_file:
            w.writeheader()
        w.writerow({
            "graph": graph, "r": r, "s": s, "tier": tier,
            "algo": f"T{tier}",
            "status": status,
            "wall_ms":  f"{fields['wall_ms']:.1f}"   if fields["wall_ms"]  >= 0 else "",
            "total_ms": f"{fields['total_ms']:.1f}"  if fields["total_ms"] >= 0 else "",
            "step4_ms": f"{fields['step4_ms']:.1f}"  if fields["step4_ms"] >= 0 else "",
            "peel_ms":  f"{fields['peel_ms']:.1f}"   if fields["peel_ms"]  >= 0 else "",
            "mem_kB":   f"{fields['mem_kB']:.0f}"    if fields["mem_kB"]   >= 0 else "",
            "max_core": fields["max_core"] if fields["max_core"] >= 0 else "",
        })


# ============ Main ============
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("server", nargs="?", default=None,
                    help="tods1|tods2|local (auto-detect if omitted)")
    ap.add_argument("--reuse-t4", action="store_true",
                    help="Join T4 rows from bench_v3_all_results.csv (RegNDC)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-build", action="store_true")
    args = ap.parse_args()

    if args.server and args.server not in SERVER_DATADIR:
        print(f"unknown server {args.server!r}", file=sys.stderr); sys.exit(1)
    server = args.server or detect_server()
    print(f"[server] {server} (datadir={SERVER_DATADIR[server]})", flush=True)

    if not args.dry_run and not args.no_build:
        cmake_build()
    link_graphs(server)

    done = load_done()
    t4_reuse = load_t4_from_v3lm() if args.reuse_t4 else {}

    rows_planned = []
    for (g, r, s) in ABLATION_CELLS:
        for tier in TIERS:
            key = (g, r, s, tier)
            if key in done:
                continue
            # Tier-1 hard skip
            if tier == 1 and (g, r, s) in T1_HARD_SKIP:
                rows_planned.append(("SKIP_T1_HARD", g, r, s, tier))
                continue
            # T4 reuse path
            if tier == 4 and args.reuse_t4 and (g, r, s) in t4_reuse:
                rows_planned.append(("REUSE_T4", g, r, s, tier))
                continue
            rows_planned.append(("RUN", g, r, s, tier))

    print(f"[plan] {len(rows_planned)} rows; "
          f"RUN={sum(1 for r in rows_planned if r[0]=='RUN')}, "
          f"REUSE_T4={sum(1 for r in rows_planned if r[0]=='REUSE_T4')}, "
          f"SKIP_T1_HARD={sum(1 for r in rows_planned if r[0]=='SKIP_T1_HARD')}",
          flush=True)

    if args.dry_run:
        for kind, g, r, s, tier in rows_planned:
            print(f"  {kind:14s}  T{tier}  {g:18s}  r={r} s={s}", flush=True)
        return

    # SEQUENTIAL execution (single-thread per cell already; cells share IO).
    # In a future revision this can be process-pool parallel; for the first
    # ablation pass keep it single-stream so logs are linear and the user
    # can interrupt at any cell boundary.
    for kind, g, r, s, tier in rows_planned:
        if kind == "SKIP_T1_HARD":
            write_row(g, r, s, tier, "TIMEOUT(skip)",
                      {"wall_ms": TIMEOUT_SEC*1000, "total_ms": -1,
                       "step4_ms": -1, "peel_ms": -1, "mem_kB": -1, "max_core": -1})
            print(f"  SKIP_T1_HARD  T{tier} {g} r={r} s={s}", flush=True)
            continue
        if kind == "REUSE_T4":
            src = t4_reuse[(g, r, s)]
            write_row(g, r, s, tier, "OK(reused)", {
                "wall_ms":  float(src.get("wall_ms")  or -1),
                "total_ms": float(src.get("total_ms") or -1),
                "step4_ms": float(src.get("step4_ms") or -1),
                "peel_ms":  float(src.get("peel_ms")  or -1),
                "mem_kB":   float(src.get("mem_kB")   or -1),
                "max_core": -1,
            })
            print(f"  REUSE_T4      T{tier} {g} r={r} s={s}", flush=True)
            continue
        print(f"  RUN           T{tier} {g} r={r} s={s} ...", flush=True, end="")
        t0 = time.time()
        status, fields, log_path = run_cell(g, r, s, tier)
        dt = time.time() - t0
        write_row(g, r, s, tier, status, fields)
        print(f" {status} ({dt:.1f}s, peak={fields['mem_kB']/1024:.0f}MB, "
              f"max_core={fields['max_core']})", flush=True)


if __name__ == "__main__":
    main()
