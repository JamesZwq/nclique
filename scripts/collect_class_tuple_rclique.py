#!/usr/bin/env python3
"""
Collect (graph, r, classes, tuples, r-cliques) for fixed s.
Append to CSV in experiments/vsafe_phase_b1/class_tuple_rclique.csv.

Reads the framework's stdout output for these counters:
  Overlap classes: N
  Active r-tuples: N
  Active r-cliques: N
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path

BIN = "./build/bin/degeneracy_cliques"
DEFAULT_OUT = "experiments/vsafe_phase_b1/class_tuple_rclique.csv"

PAT_CLASSES   = re.compile(r"Overlap classes:\s+(\d+)")
PAT_TUPLES    = re.compile(r"Active r-tuples:\s+(\d+)")
PAT_RCLIQUES  = re.compile(r"Active r-cliques:\s+(\d+)")


def run_one(graph: str, r: int, s: int, timeout_s: int = 120) -> dict:
    edges = f"graphs/{graph}.edges"
    if not os.path.exists(edges):
        return {"error": f"missing edges file: {edges}"}

    env = os.environ.copy()
    env["PIVOTER_RUN_REGION_V3LM"] = "1"
    env["PIVOTER_VSAFE_CLOUD"] = "1"
    cmd = [BIN, edges, str(r), str(s), "degen"]
    t0 = time.time()
    try:
        r_proc = subprocess.run(
            cmd, env=env, capture_output=True, text=True,
            timeout=timeout_s, check=False)
    except subprocess.TimeoutExpired:
        return {"error": f"timeout after {timeout_s}s"}
    dt = time.time() - t0
    out = r_proc.stdout + "\n" + r_proc.stderr

    cls = tup = rc = None
    m = PAT_CLASSES.search(out);  cls = int(m.group(1)) if m else None
    m = PAT_TUPLES.search(out);   tup = int(m.group(1)) if m else None
    m = PAT_RCLIQUES.search(out); rc  = int(m.group(1)) if m else None
    return {"classes": cls, "tuples": tup, "rcliques": rc, "wall_s": dt}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graphs", nargs="+", required=True,
                    help="graph names (without .edges)")
    ap.add_argument("--rs", type=int, nargs="+", default=[3, 4, 5, 6, 7, 8, 9])
    ap.add_argument("--s", type=int, default=10)
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--out", default=DEFAULT_OUT)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    new_file = not out.exists()
    fieldnames = ["graph", "r", "s", "classes", "tuples", "rcliques",
                  "wall_s", "error"]
    with out.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        if new_file:
            w.writeheader()

        for g in args.graphs:
            for r in args.rs:
                if r >= args.s:
                    continue
                print(f"  [{g} r={r} s={args.s}] running... ", end="", flush=True)
                res = run_one(g, r, args.s, args.timeout)
                row = {"graph": g, "r": r, "s": args.s,
                       "classes": res.get("classes"),
                       "tuples": res.get("tuples"),
                       "rcliques": res.get("rcliques"),
                       "wall_s": f"{res.get('wall_s', 0):.1f}",
                       "error": res.get("error", "")}
                w.writerow(row)
                f.flush()
                if res.get("error"):
                    print(f"ERR: {res['error']}")
                else:
                    print(f"cls={row['classes']} tup={row['tuples']} "
                          f"rc={row['rcliques']} ({row['wall_s']}s)")


if __name__ == "__main__":
    main()
