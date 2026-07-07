#!/usr/bin/env python3
# §131 BASELINE: the whole t=1 diagonal (r, r+1) for all r, as INDEPENDENT native runs of the
# current engine (no sharing, no host-1 shortcut). Serial; per-row TIMEOUT + RSS-poll KILL +
# engine-side SCT_MAX_INC pattern-explosion guard (exit 7) so a row can never OOM the machine.
# Protocol: ascend r=1,2,... until MAXTO consecutive aborts (r_lo edge), then descend from RTOP
# until MAXTO consecutive aborts (r_hi edge); the middle band is INFEASIBLE-for-current-engine.
# Output TSV: r status wall_s rss_mb regions merged_rc patterns total_rc maxcore mce enum build peel
import subprocess, sys, os, re, time, signal

BIN, GRAPH, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
TIMEOUT = float(sys.argv[4]) if len(sys.argv) > 4 else 240.0
RTOP    = int(sys.argv[5]) if len(sys.argv) > 5 else 130
RSSCAP_MB = float(sys.argv[6]) if len(sys.argv) > 6 else 12000.0   # hard kill above this
MAXINC  = os.environ.get("SCT_MAX_INC", "120000000")               # engine-side clean abort
MAXTO   = 2

def run_row(r):
    env = dict(os.environ); env["SCT_MAX_INC"] = MAXINC
    t0 = time.time()
    p = subprocess.Popen([BIN, GRAPH, str(r), str(r + 1)], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True, env=env)
    peak = 0.0; status = None
    while True:
        rc = p.poll()
        if rc is not None: break
        if time.time() - t0 > TIMEOUT:
            p.kill(); status = "TIMEOUT"; break
        try:
            rss = int(subprocess.run(["ps", "-o", "rss=", "-p", str(p.pid)],
                                     capture_output=True, text=True).stdout.strip() or 0) / 1024.0
            peak = max(peak, rss)
            if rss > RSSCAP_MB: p.kill(); status = "RSSKILL"; break
        except ValueError:
            pass
        time.sleep(1.0)
    try:
        o = p.communicate(timeout=30)[0] or ""
    except subprocess.TimeoutExpired:
        p.kill(); o = ""
    wall = time.time() - t0
    if status is None:
        rc = p.returncode
        status = "OK" if rc == 0 else ("EXPLODED" if rc == 7 else f"CRASH({rc})")
    d = {"status": status, "wall": wall, "rss_mb": peak}
    if "no region >= s" in o: d["status"] = "NOREGION"
    m = re.search(r"regions\(>=s\)=(\d+)", o);               d["regions"] = m.group(1) if m else "0"
    m = re.search(r"r-mergeable: (\d+) regions direct \((\d+) r-cliques\); active=(\d+)", o)
    d["merged_rc"] = m.group(2) if m else "0"
    m = re.search(r"patterns=(\d+)\s+r-cliques=(\d+)", o);   d["patterns"], d["total_rc"] = (m.group(1), m.group(2)) if m else ("0", "0")
    m = re.search(r"Max core: (\d+)", o);                    d["maxcore"] = m.group(1) if m else ""
    m = re.search(r"TIMING MCE=([\d.]+) enum=([\d.]+) sct-build\+maps=([\d.]+) peel=([\d.]+)", o)
    d["mce"], d["enum"], d["build"], d["peel"] = m.groups() if m else ("", "", "", "")
    return d

def emit(f, r, d):
    f.write(f"{r}\t{d['status']}\t{d['wall']:.2f}\t{d.get('rss_mb',0):.0f}\t{d.get('regions','')}\t"
            f"{d.get('merged_rc','')}\t{d.get('patterns','')}\t{d.get('total_rc','')}\t{d.get('maxcore','')}\t"
            f"{d.get('mce','')}\t{d.get('enum','')}\t{d.get('build','')}\t{d.get('peel','')}\n")
    f.flush()
    print(f"r={r}: {d['status']} {d['wall']:.1f}s rss={d.get('rss_mb',0):.0f}MB pats={d.get('patterns','')}", flush=True)

ABORTS = ("TIMEOUT", "RSSKILL", "EXPLODED", "CRASH")
with open(OUT, "w") as f:
    f.write("r\tstatus\twall_s\trss_mb\tregions\tmerged_rc\tpatterns\ttotal_rc\tmaxcore\tmce\tenum\tbuild\tpeel\n")
    r, consec, r_lo = 1, 0, None
    while True:
        d = run_row(r); emit(f, r, d)
        if d["status"] == "NOREGION": r_lo = None; break
        if d["status"].startswith(ABORTS):
            consec += 1
            if consec >= MAXTO: r_lo = r - MAXTO + 1; break
        else: consec = 0
        r += 1
    if r_lo is not None:
        consec, r_hi = 0, None
        for rr in range(RTOP, r_lo, -1):
            d = run_row(rr)
            if d["status"] == "NOREGION": continue
            emit(f, rr, d)
            if d["status"].startswith(ABORTS):
                consec += 1
                if consec >= MAXTO: r_hi = rr + MAXTO - 1; break
            else: consec = 0
        f.write(f"# INFEASIBLE-BAND (current engine, timeout {TIMEOUT:.0f}s, rss-cap {RSSCAP_MB:.0f}MB, inc-cap {MAXINC}): r in [{r_lo}, {r_hi}]\n")
        print(f"INFEASIBLE-BAND: [{r_lo}, {r_hi}]", flush=True)
print("DIAG_BASELINE_DONE", flush=True)
