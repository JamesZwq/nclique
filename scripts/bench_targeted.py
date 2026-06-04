#!/usr/bin/env python3
"""Focused bench driver: run one algorithm on a (graph × r × s) grid.

Used for filling in specific coverage gaps in the SIGMOD 2027 paper:
  - V3LM_NOCPI ablation row in Fig 3
  - REF + V3LM at r=10, r=15 for memory ratio fig
  - Large-graph extension (com-lj, com-orkut) for Fig 3 coverage

Output: one CSV row per (graph, r, s, algo) cell with
  status / wall_ms / total_ms / step4_ms / peel_ms / hier_ms / mem_kB.

Parses the same timing markers as bench_v3_all.py.
"""
from __future__ import annotations
import argparse, csv, os, re, signal, subprocess, sys, time
from pathlib import Path

BIN = "./build/bin/degeneracy_cliques"
ALGO_ENV = {
    "REF":         "PIVOTER_RUN_REF",
    "RegNDC":      "PIVOTER_RUN_REGION_V3LM",
    "V3LM":        "PIVOTER_RUN_REGION_V3LM",
    "V3LM_NOCPI":  "PIVOTER_RUN_REGION_V3LM_NOCPI",
    "V3LM_HIER":   "PIVOTER_RUN_REGION_V3LM_HIER",
}
RE_TOTAL = re.compile(r'Total time:\s*([\d.]+)\s*ms')
RE_STEP4 = re.compile(r'(?:CPI counting time|NoCPI enumeration time):\s*([\d.]+)')
RE_PEEL  = re.compile(r'Peeling time:\s*([\d.]+)')
RE_HIER  = re.compile(r'Hierarchy time:\s*([\d.]+)')


def parse_time_v(stderr: str) -> int:
    """Extract Maximum resident set size (kB) from /usr/bin/time -v output."""
    m = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', stderr)
    return int(m.group(1)) if m else -1


def run_cell(graph: str, r: int, s: int, algo: str,
             timeout_s: int, mem_cap_kb: int) -> dict:
    edges = f"graphs/{graph}.edges"
    if not os.path.exists(edges):
        return {"status": "MISSING_EDGES"}

    env = os.environ.copy()
    env[ALGO_ENV[algo]] = "1"
    if algo in ("RegNDC", "V3LM"):
        env["PIVOTER_VSAFE_CLOUD"] = "1"

    cmd = ["/usr/bin/time", "-v", BIN, edges, str(r), str(s), "degen"]
    t0 = time.time()
    proc = subprocess.Popen(
        cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, preexec_fn=os.setsid)
    try:
        stdout, stderr = proc.communicate(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        # Kill the entire process group (covers /usr/bin/time wrapper AND binary).
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()
        return {"status": "TIMEOUT", "wall_ms": timeout_s * 1000}
    # Synthesize CompletedProcess-like view for the rest of the function.
    class _R: pass
    proc_view = _R()
    proc_view.returncode = proc.returncode
    proc_view.stdout = stdout
    proc_view.stderr = stderr
    proc = proc_view
    dt_ms = (time.time() - t0) * 1000

    stdout = proc.stdout
    stderr = proc.stderr
    mem_kb = parse_time_v(stderr)
    if mem_kb > mem_cap_kb:
        return {"status": "OOM", "wall_ms": dt_ms, "mem_kB": mem_kb}
    if proc.returncode != 0:
        # detect killed by OOM
        if "Killed" in stderr or proc.returncode == -9:
            return {"status": "OOM", "wall_ms": dt_ms, "mem_kB": mem_kb}
        return {"status": f"ERROR({proc.returncode})",
                "wall_ms": dt_ms, "mem_kB": mem_kb}

    total_ms = float(m.group(1)) if (m := RE_TOTAL.search(stdout)) else dt_ms
    step4_ms = float(m.group(1)) if (m := RE_STEP4.search(stdout)) else -1
    peel_ms  = float(m.group(1)) if (m := RE_PEEL.search(stdout))  else -1
    hier_ms  = float(m.group(1)) if (m := RE_HIER.search(stdout))  else -1

    return {
        "status": "OK",
        "wall_ms": dt_ms, "total_ms": total_ms,
        "step4_ms": step4_ms, "peel_ms": peel_ms, "hier_ms": hier_ms,
        "mem_kB": mem_kb,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", nargs="+", required=True)
    ap.add_argument("--rs", type=int, nargs="+", required=True)
    ap.add_argument("--ss", type=int, nargs="+", required=True)
    ap.add_argument("--algo", required=True, choices=list(ALGO_ENV))
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--mem-cap-gb", type=int, default=250)
    ap.add_argument("--out", required=True, help="path to append CSV rows")
    ap.add_argument("--skip-after-fail", action="store_true",
                    help="for fixed r, skip higher s once a cell fails")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not out_path.exists()
    fieldnames = ["graph","r","s","algo","status",
                  "wall_ms","total_ms","step4_ms","peel_ms","hier_ms","mem_kB"]

    # Load existing cells to skip: any non-blank status counts as already-tried,
    # so a restart does not retry TIMEOUT / OOM / SKIP_AFTER_FAIL cells.
    done = set()
    if out_path.exists():
        with out_path.open() as f:
            for row in csv.DictReader(f):
                if row.get("status", "").strip():
                    done.add((row["graph"], int(row["r"]),
                              int(row["s"]), row["algo"]))

    mem_cap_kb = args.mem_cap_gb * 1024 * 1024

    with out_path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if new_file:
            w.writeheader()
            f.flush()

        for g in args.graphs:
            for r in sorted(args.rs):
                failed_at_s = None
                for s in sorted(args.ss):
                    if r >= s:
                        continue
                    if args.skip_after_fail and failed_at_s is not None and s > failed_at_s:
                        w.writerow({"graph":g,"r":r,"s":s,"algo":args.algo,
                                    "status":"SKIP_AFTER_FAIL","wall_ms":"",
                                    "total_ms":"","step4_ms":"","peel_ms":"",
                                    "hier_ms":"","mem_kB":""})
                        f.flush()
                        continue
                    key = (g, r, s, args.algo)
                    if key in done:
                        continue
                    print(f"  [{g} r={r} s={s} {args.algo}] running... ",
                          end="", flush=True)
                    res = run_cell(g, r, s, args.algo, args.timeout, mem_cap_kb)
                    row = {"graph":g,"r":r,"s":s,"algo":args.algo,
                           **{k: f"{v:.1f}" if isinstance(v,float) else str(v)
                              for k,v in res.items()}}
                    for k in fieldnames:
                        row.setdefault(k, "")
                    w.writerow(row)
                    f.flush()
                    status = res.get("status","?")
                    wall = res.get("wall_ms",0)
                    mem = res.get("mem_kB",0)
                    print(f"{status}  {wall:.0f}ms  {mem/1024:.0f}MB")
                    if status in ("TIMEOUT","OOM","MISSING_EDGES"):
                        failed_at_s = s
    print(f"\ndone. CSV at {out_path}")


if __name__ == "__main__":
    main()
