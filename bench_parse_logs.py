#!/usr/bin/env python3
"""Parse bench_v3_phases log files → clean CSV + summary tables."""

import os, re, csv, sys
from pathlib import Path
from collections import defaultdict

LOGDIR = Path("bench_v3_phases_logs")
OUTCSV = Path("bench_v3_phases_results.csv")

# Regex patterns for extracting metrics from log output
PATTERNS = {
    # Timing
    "total_ms":      r'NucleusCoreDecomposition took:\s*([\d.]+)\s*ms',
    "sdct_ms":       r'SDCT_(?:MaxClique|Fused) took:\s*([\d.]+)\s*ms',
    "maxcliq_ms":    r'MaxCliqEnum \(V3\) took:\s*([\d.]+)\s*ms',
    "maxcliq_v2_ms": r'MaxCliqEnum took:\s*([\d.]+)\s*ms',
    "ci_build_ms":   r'cliqueIndex mapList build took:\s*([\d.]+)\s*ms',
    "merge_ms":      r'r-Mergeable classification:\s*([\d.]+)\s*ms',
    "cpi_ms":        r'CPI counting time:\s*([\d.]+)',
    "pathinfo_ms":   r'PathInfo build time:\s*([\d.]+)',
    "peel_ms":       r'Peeling time:\s*([\d.]+)',
    "v3_total_ms":   r'Total time:\s*([\d.]+)',
    "st_fused_ms":   r'fused build\+counting \(ST\) took:\s*([\d.]+)',

    # Memory (kB)
    "graph_mem_kB":  r'\[Memory-\w+\]\s*Graph Memory:\s*([\d.]+)\s*kB',
    "index_mem_kB":  r'\[Memory-\w+\]\s*Other Index Memory:\s*([\d.]+)\s*kB',
    "final_mem_kB":  r'\[Memory-\w+\]\s*Final Memory:\s*([\d.]+)\s*kB',
    "rclique_mem_kB":r'\[Memory-\w+\]\s*r-clique index:\s*([\d.]+)\s*kB',

    # Stats
    "leaf_count":    r'SDCT leaf count(?: \(stored\))?:\s*(\d+)',
    "max_core":      r'Max core:\s*(\d+)',
    "rcliques":      r'r-cliques:\s*(\d+)',
    "fully_merge":   r'Fully mergeable regions:\s*(\d+)',
    "remaining":     r'Remaining regions:\s*(\d+)',
    "classes":       r'Overlap classes:\s*(\d+)',
    "tuples":        r'r-tuples:\s*(\d+)',
    "recursive":     r'Total recursive calls:\s*(\d+)',
    "v2_stuples":    r's-tuples:\s*(\d+)',
    "maxcliq_count": r'maximal cliques \((?:minSize|≥s)=\d+\):\s*(\d+)',
    "maxcliq_count2":r'MaxCliqEnum:\s*(\d+)\s*maximal cliques',

    # Correctness
    "cpi_match":     r'CPI vs CPath match:\s*(\w+)',
    "pathinfo_match":r'CPI vs PathInfo match:\s*(\w+)',
}

FIELDS = [
    "graph", "r", "s", "algo", "status",
    "total_ms", "sdct_ms", "maxcliq_ms", "maxcliq_v2_ms", "ci_build_ms",
    "merge_ms", "cpi_ms", "pathinfo_ms", "peel_ms", "v3_total_ms", "st_fused_ms",
    "graph_mem_kB", "index_mem_kB", "final_mem_kB", "rclique_mem_kB",
    "leaf_count", "max_core", "rcliques",
    "fully_merge", "remaining", "classes", "tuples", "recursive",
    "v2_stuples", "maxcliq_count", "maxcliq_count2",
    "cpi_match", "pathinfo_match",
]


def parse_log(logfile: Path) -> dict:
    """Parse a single log file and extract all metrics."""
    text = logfile.read_text(errors="replace")
    row = {}
    for key, pattern in PATTERNS.items():
        m = re.search(pattern, text)
        row[key] = m.group(1) if m else ""
    return row


def parse_filename(name: str):
    """Extract graph, r, s, algo from filename like 'com-dblp_r3_s4_ST.log'."""
    m = re.match(r'^(.+)_r(\d+)_s(\d+)_(\w+)\.log$', name)
    if not m:
        return None
    return m.group(1), m.group(2), m.group(3), m.group(4)


def detect_status(logfile: Path, row: dict) -> str:
    """Determine job status from log content."""
    text = logfile.read_text(errors="replace")
    if row.get("total_ms"):
        return "OK"
    if "OOM" in text or "Cannot allocate memory" in text or "bad_alloc" in text:
        return "OOM"
    if "Killed" in text or "signal 9" in text:
        return "KILLED"
    if "Segmentation fault" in text or "signal 11" in text:
        return "SEGFAULT"
    if "TIMEOUT" in text:
        return "TIMEOUT"
    # Check if file is very small (job was killed before output)
    if logfile.stat().st_size < 100:
        return "KILLED"
    return "ERROR"


def build_csv():
    """Scan all log files and build CSV."""
    if not LOGDIR.exists():
        print(f"Error: {LOGDIR} not found")
        sys.exit(1)

    rows = []
    for logfile in sorted(LOGDIR.glob("*.log")):
        parsed = parse_filename(logfile.name)
        if not parsed:
            continue
        graph, r, s, algo = parsed

        row = parse_log(logfile)
        row["graph"] = graph
        row["r"] = r
        row["s"] = s
        row["algo"] = algo
        row["status"] = detect_status(logfile, row)
        rows.append(row)

    # Write CSV
    with open(OUTCSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {len(rows)} rows to {OUTCSV}")
    return rows


def print_summary(rows):
    """Print summary tables."""
    # Group by (graph, r, s)
    data = defaultdict(dict)
    for row in rows:
        key = (row["graph"], int(row["r"]), int(row["s"]))
        algo = row["algo"]
        status = row["status"]
        total = row.get("total_ms", "")
        mem = row.get("final_mem_kB", "")
        if status == "OK" and total:
            data[key][algo] = {"time": float(total), "mem": float(mem) if mem else None, "status": "OK"}
        else:
            data[key][algo] = {"time": None, "mem": None, "status": status}

    # Sort keys
    keys = sorted(data.keys())

    # Time table
    print("\n" + "=" * 80)
    print("  Total Time (seconds)")
    print("=" * 80)
    print(f"{'Graph':<18} {'r':>2} {'s':>2} | {'ST':>10} {'V2':>10} {'V3':>10} | Winner")
    print("-" * 80)

    for key in keys:
        graph, r, s = key
        d = data[key]

        def fmt(algo):
            info = d.get(algo, {"time": None, "status": "—"})
            if info["time"] is not None:
                return f"{info['time']/1000:.1f}s"
            return info["status"]

        vals = {}
        for algo in ["ST", "V2", "V3"]:
            info = d.get(algo, {"time": None, "status": "—"})
            vals[algo] = info

        # Find winner
        times = {a: v["time"] for a, v in vals.items() if v["time"] is not None}
        if times:
            winner = min(times, key=times.get)
            # Calculate speedup vs second best
            sorted_times = sorted(times.values())
            if len(sorted_times) >= 2:
                speedup = sorted_times[1] / sorted_times[0]
                winner_str = f"{winner} ({speedup:.1f}x)"
            else:
                winner_str = winner
        else:
            winner_str = "—"

        print(f"{graph:<18} {r:>2} {s:>2} | {fmt('ST'):>10} {fmt('V2'):>10} {fmt('V3'):>10} | {winner_str}")

    # Memory table
    print("\n" + "=" * 80)
    print("  Peak Memory (MB)")
    print("=" * 80)
    print(f"{'Graph':<18} {'r':>2} {'s':>2} | {'ST':>10} {'V2':>10} {'V3':>10}")
    print("-" * 80)

    for key in keys:
        graph, r, s = key
        d = data[key]

        def fmt_mem(algo):
            info = d.get(algo, {"mem": None, "status": "—"})
            if info["mem"] is not None:
                return f"{info['mem']/1024:.0f}MB"
            return info["status"]

        print(f"{graph:<18} {r:>2} {s:>2} | {fmt_mem('ST'):>10} {fmt_mem('V2'):>10} {fmt_mem('V3'):>10}")


def main():
    rows = build_csv()
    if rows:
        print_summary(rows)


if __name__ == "__main__":
    main()
