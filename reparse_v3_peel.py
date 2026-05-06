#!/usr/bin/env python3
"""
Re-parse build_ms and peel_ms from existing bench_r1_v3 logs into the CSV.

The first cut of bench_r1_v3.py only stored took_ms = build + peel as a
single column. paper §7.2 quotes peel-only speedup, so we need build_ms
and peel_ms as separate columns. Logs are full stdout+stderr per cell —
re-parse them and add the new columns to the existing CSV.

Idempotent: rows that already have build_ms/peel_ms are left alone.

Usage:
    python3 reparse_v3_peel.py
    # rewrites paper_data/01_main_benchmark_v3.csv in place,
    # backing up the original to .csv.bak first.
"""
from __future__ import annotations
import csv, re, shutil, sys
from pathlib import Path

CSV_PATH = Path("paper_data/01_main_benchmark_v3.csv")
LOG_DIR  = Path("bench_r1_v3_logs")
BACKUP   = CSV_PATH.with_suffix(".csv.bak")

_LOAD_RE  = re.compile(r"loadAndSort took:\s*([\d.]+)\s*ms")
_BUILD_RE = re.compile(r"(?:buildSDCT|SDCT_Fused|SDCT\+callback) took:\s*([\d.]+)\s*ms")
_PEEL_V3_RE = re.compile(r"ST_V3 r=1 \(peel\) took:\s*([\d.]+)\s*ms")
_PEEL_REF_RE = re.compile(r"NucleusCoreDecomposition took:\s*([\d.]+)\s*ms")


def parse_log(path: Path, algo: str) -> tuple[str, str]:
    """Return (build_ms, peel_ms) as strings (or '') extracted from a log."""
    if not path.exists():
        return ("", "")
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return ("", "")
    build = ""; peel = ""
    if (m := _BUILD_RE.search(text)): build = m.group(1)
    if algo == "Pure":
        if (m := _PEEL_V3_RE.search(text)): peel = m.group(1)
        elif (m := _PEEL_REF_RE.search(text)): peel = m.group(1)
    else:
        if (m := _PEEL_REF_RE.search(text)): peel = m.group(1)
    return (build, peel)


def main():
    if not CSV_PATH.exists():
        sys.exit(f"missing {CSV_PATH}")
    rows = list(csv.DictReader(CSV_PATH.open()))
    if not rows:
        sys.exit("no rows in CSV")
    fieldnames = list(rows[0].keys())
    # Insert build_ms / peel_ms after took_ms if not already present.
    if "build_ms" not in fieldnames:
        i = fieldnames.index("took_ms") + 1
        fieldnames.insert(i, "build_ms")
        fieldnames.insert(i+1, "peel_ms")

    n_filled = 0
    n_already = 0
    n_missing_log = 0
    for r in rows:
        if r.get("build_ms") or r.get("peel_ms"):
            n_already += 1
            continue
        if r.get("status") != "OK":
            continue
        graph = r.get("graph", "")
        s     = r.get("s", "")
        algo  = r.get("algorithm", "")
        run   = r.get("run", "0")
        log = LOG_DIR / f"{graph}_s{s}_{algo}_r{run}.log"
        b, p = parse_log(log, algo)
        if not (b or p):
            n_missing_log += 1
            continue
        r["build_ms"] = b
        r["peel_ms"]  = p
        n_filled += 1

    # Backup once.
    if not BACKUP.exists():
        shutil.copy2(CSV_PATH, BACKUP)
        print(f"[backup] {CSV_PATH} -> {BACKUP}")
    with CSV_PATH.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            # Make sure all rows have all keys (DictWriter is strict).
            for k in fieldnames: r.setdefault(k, "")
            w.writerow(r)
    print(f"[reparse] filled={n_filled}  already_had={n_already}  no_log={n_missing_log}")


if __name__ == "__main__":
    main()
