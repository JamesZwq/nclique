#!/usr/bin/env python3
"""
bench_region_native.py — server comparison of the region-native counting
engine (region_native/region_native.cpp) vs the existing CPI/SDCT engine,
for the "fuse pivot into Region/Class" milestone (SigmodPlus 24.2).

Measures, per (graph, r, s), the cost of GETTING TO per-tuple support:
  region-native: degeneracy MCE + region union-count support
  existing CPI : SDCT build + MaxCliqEnum + CPI counting + PathInfo build
Both produce the same per-tuple support (bit-exact, validated offline vs
direct enumeration). This isolates whether dropping the pivoted SDCT in
favour of a region-native count is faster; the headline is web-it, where
the SDCT build dominates.

Also runs region_native --verify with a small sample per cell as an
on-server correctness gate.

Usage: nohup python3 bench_region_native.py > brn.log 2>&1 &
"""
from __future__ import annotations
import csv, os, re, subprocess, time
from pathlib import Path

RN  = "./region_native/region_native"
BIN = "./build/bin/degeneracy_cliques"
GRAPHS = ["dblp-core30", "ca-GrQc", "ca-HepPh", "ca-CondMat", "com-dblp", "web-it-2004"]
# FULL (r,s) grid to test the size-free criterion: region-native time must
# NOT grow significantly with r or s. region-native runs every cell (fast);
# CPI uses a skip-floor (once it times out at some s, higher s same r is
# skipped) since it explodes with s -- that explosion vs our flatness IS
# the result.
RGRID = [int(x) for x in os.getenv("BRN_RGRID", "3,4,5,6,7").split(",")]
SMAX = int(os.getenv("BRN_SMAX", "20"))
RS = [(r, s) for r in RGRID for s in range(r + 1, SMAX + 1)]
TIMEOUT = int(os.getenv("BRN_TIMEOUT", "1800"))
VERIFY = int(os.getenv("BRN_VERIFY", "2000"))
OUT = Path(os.getenv("BRN_OUT", "paper_data/bench_region_native.csv"))
LOGD = Path("brn_logs"); LOGD.mkdir(exist_ok=True)
FIELDS = ["graph", "r", "s", "regions", "tuples",
          "rn_mce_s", "rn_support_s", "rn_total_s", "rn_verify",
          "cpi_sdct_s", "cpi_mcenum_s", "cpi_count_s", "cpi_pathinfo_s",
          "cpi_total_s", "speedup"]

def grep1(txt, pat, grp=1, cast=float, default=-1):
    m = re.search(pat, txt)
    return cast(m.group(grp)) if m else default

def run(cmd, to):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=to)
        return p.stdout + p.stderr
    except subprocess.TimeoutExpired:
        return "__TIMEOUT__"

def load_done(p):
    d = set()
    if p.exists():
        for row in csv.DictReader(open(p)):
            d.add((row["graph"], int(row["r"]), int(row["s"])))
    return d

def main():
    OUT.parent.mkdir(parents=True, exist_ok=True)
    done = load_done(OUT)
    new = not OUT.exists()
    f = open(OUT, "a", newline=""); w = csv.DictWriter(f, fieldnames=FIELDS)
    if new: w.writeheader(); f.flush()
    cpi_skip = set()   # (graph, r) where CPI timed out -> skip higher s
    for g in GRAPHS:
        gp = f"graphs/{g}.edges"
        if not os.path.exists(gp): print(f"[skip] {g}", flush=True); continue
        for (r, s) in RS:
            if (g, r, s) in done: continue
            # region-native (with small verify sample)
            rn = run([RN, gp, str(r), str(s), "--verify", str(VERIFY),
                      "--mce-budget", str(TIMEOUT)], TIMEOUT + 60)
            (LOGD / f"{g}_{r}_{s}_rn.log").write_text(rn if rn != "__TIMEOUT__" else "TIMEOUT")
            if rn == "__TIMEOUT__":
                row = dict(graph=g, r=r, s=s, regions=-1, tuples=-1, rn_mce_s=-1,
                           rn_support_s=-1, rn_total_s=-1, rn_verify="TIMEOUT",
                           cpi_sdct_s=-1, cpi_mcenum_s=-1, cpi_count_s=-1,
                           cpi_pathinfo_s=-1, cpi_total_s=-1, speedup=-1)
                w.writerow(row); f.flush(); print(row, flush=True); continue
            regions = grep1(rn, r"regions\(>=s\)=(\d+)", cast=int)
            tuples = grep1(rn, r"tuples=(\d+)", cast=int)
            rn_mce = grep1(rn, r"MCE=([\d.]+)s support")
            rn_sup = grep1(rn, r"support=([\d.]+)s \(support-only")
            verify = "EXACT" if "[EXACT]" in rn else ("MISMATCH" if "MISMATCH" in rn else "?")
            rn_total = (rn_mce + rn_sup) if (rn_mce >= 0 and rn_sup >= 0) else -1
            # existing CPI engine (skip-floor: CPI explodes with s)
            if (g, r) in cpi_skip:
                cp = "__SKIP__"
            else:
                os.environ["PIVOTER_RUN_REGION_V3LM"] = "1"
                cp = run([BIN, gp, str(r), str(s), "degen"], TIMEOUT + 60)
                if cp == "__TIMEOUT__":
                    cpi_skip.add((g, r))
            (LOGD / f"{g}_{r}_{s}_cpi.log").write_text(cp if cp != "__TIMEOUT__" else "TIMEOUT")
            sdct = grep1(cp, r"SDCT_(?:MaxClique|Fused) took: ([\d.]+) ms") / 1000.0
            mcenum = grep1(cp, r"MaxCliqEnum \(V3/V4\) took: ([\d.]+) ms") / 1000.0
            ccount = grep1(cp, r"CPI counting time: ([\d.]+) ms") / 1000.0
            pinfo = grep1(cp, r"PathInfo build time: ([\d.]+) ms") / 1000.0
            cpi_total = sum(x for x in (sdct, mcenum, ccount, pinfo) if x >= 0) \
                        if cp not in ("__TIMEOUT__","__SKIP__") and sdct >= 0 else -1
            speedup = (cpi_total / rn_total) if (cpi_total > 0 and rn_total > 0) else -1
            row = dict(graph=g, r=r, s=s, regions=regions, tuples=tuples,
                       rn_mce_s=round(rn_mce, 3), rn_support_s=round(rn_sup, 3),
                       rn_total_s=round(rn_total, 3), rn_verify=verify,
                       cpi_sdct_s=round(sdct, 3), cpi_mcenum_s=round(mcenum, 3),
                       cpi_count_s=round(ccount, 3), cpi_pathinfo_s=round(pinfo, 3),
                       cpi_total_s=round(cpi_total, 3) if cpi_total >= 0 else -1,
                       speedup=round(speedup, 2) if speedup > 0 else -1)
            w.writerow(row); f.flush()
            print(f"{g:13} ({r},{s}) RN={rn_total:.2f}s CPI={cpi_total:.2f}s "
                  f"speedup={speedup:.2f}x verify={verify}", flush=True)
    f.close()
    print("[brn] done", flush=True)

if __name__ == "__main__":
    main()
