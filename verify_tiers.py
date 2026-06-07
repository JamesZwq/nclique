#!/usr/bin/env python3
"""
verify_tiers.py — smoke verification that T1/T2/T3/T4 produce the same
per-tuple κ_s histogram on a small graph.

The κ_s histogram is the multiset of (core, count) pairs printed by the
binary as the "Core value distribution" block.  All four tiers should
produce identical histograms when run on the same (graph, r, s).

Usage:
  python3 verify_tiers.py            # default: ca-GrQc at r=3 s=4
  python3 verify_tiers.py ca-GrQc 5 6
  python3 verify_tiers.py --all      # run a small suite (paper-example,
                                     # bio-celegans, ca-GrQc on 3 cells)
"""
from __future__ import annotations
import argparse, os, re, subprocess, sys
from pathlib import Path

BIN = Path("./build/bin/degeneracy_cliques")
TIMEOUT = 600


def run_tier(graph: str, r: int, s: int, tier: int) -> dict:
    env = os.environ.copy()
    env["PIVOTER_RUN_REGION_TIER"] = "1"
    env["PIVOTER_TIER"] = str(tier)
    cmd = [str(BIN), f"graphs/{graph}.edges", str(r), str(s), "degen"]
    try:
        p = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {"status": "TIMEOUT"}
    txt = p.stdout + p.stderr
    # Parse `core=N count=M` histogram lines printed by the binary.
    dist = {}
    for line in txt.splitlines():
        m = re.match(r'\s*core=(\d+(?:\.\d+)?)\s+count=(\d+)', line)
        if m:
            core = int(float(m.group(1)))
            cnt  = int(m.group(2))
            dist[core] = dist.get(core, 0) + cnt
    m_max = re.search(r'Max core:\s*(\d+)', txt)
    m_total = re.search(r'Total time:\s*(\d+)', txt)
    return {
        "status": "OK" if p.returncode == 0 else f"ERR({p.returncode})",
        "max_core": int(m_max.group(1)) if m_max else -1,
        "total_ms": int(m_total.group(1)) if m_total else -1,
        "dist": dist,
        "dist_sum": sum(dist.values()),
        "stdout_path": None,
    }


def compare_cell(graph: str, r: int, s: int) -> bool:
    print(f"\n=== {graph} r={r} s={s} ===", flush=True)
    results = {}
    for tier in (1, 2, 3, 4):
        results[tier] = run_tier(graph, r, s, tier)
        rr = results[tier]
        print(f"  T{tier}: status={rr['status']:8s}  "
              f"max_core={rr['max_core']:4d}  "
              f"dist_keys={len(rr['dist'])}  "
              f"dist_sum={rr['dist_sum']}  "
              f"total={rr['total_ms']}ms",
              flush=True)
    ok = True
    ref = results[4]
    if ref["status"] != "OK":
        print(f"  ! T4 not OK ({ref['status']}); cannot compare", flush=True)
        return False
    for tier in (1, 2, 3):
        rr = results[tier]
        if rr["status"] != "OK":
            print(f"  ! T{tier} not OK ({rr['status']})", flush=True)
            ok = False; continue
        if rr["max_core"] != ref["max_core"]:
            print(f"  ! T{tier} max_core={rr['max_core']} != T4={ref['max_core']}",
                  flush=True)
            ok = False
        if rr["dist"] != ref["dist"]:
            # Detailed diff
            diff_keys = set(rr["dist"]) ^ set(ref["dist"])
            mismatch = [(k, rr["dist"].get(k, 0), ref["dist"].get(k, 0))
                        for k in sorted(set(rr["dist"]) | set(ref["dist"]))
                        if rr["dist"].get(k, 0) != ref["dist"].get(k, 0)]
            print(f"  ! T{tier} dist != T4 dist; first 5 mismatches: {mismatch[:5]}",
                  flush=True)
            ok = False
    print(f"  -> {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("graph", nargs="?", default="ca-GrQc")
    ap.add_argument("r", nargs="?", type=int, default=3)
    ap.add_argument("s", nargs="?", type=int, default=4)
    ap.add_argument("--all", action="store_true",
                    help="Run a small suite across multiple cells")
    args = ap.parse_args()

    cells = []
    if args.all:
        for g, r, s in [
            ("bio-celegans", 3, 4),
            ("ca-GrQc",      3, 4),
            ("ca-GrQc",      4, 5),
            ("ca-GrQc",      5, 6),
            ("ca-CondMat",   3, 4),
        ]:
            if Path(f"graphs/{g}.edges").exists():
                cells.append((g, r, s))
    else:
        cells.append((args.graph, args.r, args.s))

    if not cells:
        print("No graphs found in ./graphs/", file=sys.stderr); sys.exit(1)

    all_pass = True
    for (g, r, s) in cells:
        if not compare_cell(g, r, s):
            all_pass = False

    print("\n" + ("ALL PASS" if all_pass else "SOME FAILED"))
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
