#!/usr/bin/env python3
"""
CS10 hierarchy via paper Def 4 (s-connected nuclei).

Replaces build_hierarchy.py (original used edge connectivity, which only
matches Def 4 at s=2).  Drives the V3 binary with PIVOTER_DUMP_HIER and
collects:

    16_cs10_hierarchy_summary.csv   per-(graph, s) aggregate stats
    17_cs10_persistence.csv         per-hierarchy-node row

The binary's hier dump (id,k_birth,k_death,parent,size_birth,size_death,
persistence) is post-processed here into these two paper CSVs.

Usage:
    python3 build_hierarchy_v3.py
"""
from __future__ import annotations
import csv, os, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BIN  = ROOT / "build/bin/degeneracy_cliques"

GRAPHS = [
    ("com-dblp",      ROOT / "graphs/com-dblp.edges"),
    ("com-youtube",   ROOT / "graphs/com-youtube.edges"),
    ("web-Stanford",  ROOT / "graphs/web-Stanford.edges"),
]
S_VALUES = [3, 5, 10, 15]
MIN_SIZE = 2

OUT_SUM  = ROOT / "paper_data/16_cs10_hierarchy_summary.csv"
OUT_PERS = ROOT / "paper_data/17_cs10_persistence.csv"


def run_one(graph_path: Path, s: int) -> tuple[list[dict], int]:
    """Run V3 hier dump.  Returns (hier_rows, nz_vertices)."""
    with tempfile.TemporaryDirectory() as td:
        hier = Path(td) / "hier.csv"
        core = Path(td) / "core.tsv"
        env = os.environ.copy()
        env["PIVOTER_RUN_ST_V3"]    = "1"
        env["PIVOTER_DUMP_HIER"]    = str(hier)
        env["PIVOTER_DUMP_CORE"]    = str(core)
        env["PIVOTER_HIER_MIN_SIZE"] = str(MIN_SIZE)
        env["OMP_NUM_THREADS"]      = "1"
        proc = subprocess.run(
            [str(BIN), str(graph_path), "1", str(s)],
            env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stderr)
            sys.exit(f"binary failed on {graph_path.name} s={s}")
        rows = list(csv.DictReader(hier.open()))
        # Count vertices with core > 0; skip '#' header lines.
        nz = 0
        with core.open() as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        if float(parts[1]) > 0: nz += 1
                    except ValueError:
                        pass
        return rows, nz


def summarize(rows: list[dict]) -> dict:
    """Compute aggregate stats for one (graph, s) cell."""
    if not rows:
        return dict(n_nodes=0, n_leaves=0, n_internal=0, max_depth=0,
                    max_persistence=0, total_persistence=0)
    children = defaultdict(list)
    by_id = {}
    for r in rows:
        nid = int(r["id"])
        parent = int(r["parent"])
        by_id[nid] = r
        if parent != -1:
            children[parent].append(nid)
    n_nodes = len(rows)
    n_leaves = sum(1 for nid in by_id if not children[nid])
    n_internal = n_nodes - n_leaves
    # iterative depth via post-order traversal
    depth = {}
    order = []
    visited = set()
    for root in [int(r["id"]) for r in rows if int(r["parent"]) == -1]:
        stack = [(root, False)]
        while stack:
            n, done = stack.pop()
            if done:
                depth[n] = 1 + max((depth[c] for c in children[n]), default=0)
                order.append(n)
                continue
            if n in visited: continue
            visited.add(n)
            stack.append((n, True))
            for c in children[n]:
                stack.append((c, False))
    max_depth = max(depth.values(), default=0)
    persistences = [int(r["persistence"]) for r in rows]
    return dict(n_nodes=n_nodes, n_leaves=n_leaves, n_internal=n_internal,
                max_depth=max_depth,
                max_persistence=max(persistences, default=0),
                total_persistence=sum(persistences))


def main():
    if not BIN.exists():
        sys.exit(f"binary not found: {BIN}.  cmake --build build first")
    OUT_SUM.parent.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    persistence_rows = []
    for gname, gpath in GRAPHS:
        for s in S_VALUES:
            print(f"[run] {gname} s={s}", flush=True)
            rows, nz = run_one(gpath, s)
            stats = summarize(rows)
            summary_rows.append({
                "graph": gname, "s": s, "nz_vertices": nz,
                **stats,
            })
            for r in rows:
                persistence_rows.append({
                    "graph": gname, "s": s,
                    "node_id":      r["id"],
                    "k_birth":      r["k_birth"],
                    "k_death":      r["k_death"],
                    "persistence":  r["persistence"],
                    "size_birth":   r["size_birth"],
                    "size_death":   r["size_death"],
                })
            print(f"  -> nz={nz} nodes={stats['n_nodes']} "
                  f"max_depth={stats['max_depth']} max_pers={stats['max_persistence']}",
                  flush=True)

    # Write summary
    with OUT_SUM.open("w", newline="") as f:
        cols = ["graph", "s", "nz_vertices", "n_nodes", "n_leaves",
                "n_internal", "max_depth", "max_persistence", "total_persistence"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in summary_rows: w.writerow(r)
    print(f"wrote {OUT_SUM}")

    # Write persistence
    with OUT_PERS.open("w", newline="") as f:
        cols = ["graph", "s", "node_id", "k_birth", "k_death",
                "persistence", "size_birth", "size_death"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in persistence_rows: w.writerow(r)
    print(f"wrote {OUT_PERS}  ({len(persistence_rows)} rows)")


if __name__ == "__main__":
    main()
