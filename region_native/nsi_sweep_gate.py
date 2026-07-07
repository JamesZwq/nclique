#!/usr/bin/env python3
# §129 NSI/FPS sweep gate: run the sweep binary once (SCT_SWEEP=smax) and the native binary once
# per cell; compare per-cell core distributions EXACTLY. The sweep's core=0 line is dropped first:
# it counts r-cliques in no s-clique (hosted only in regions < s), which the native per-cell run
# never enumerates (its MCE floor is s). Everything else must match line-for-line.
import subprocess, sys, os, re, time

TIMEV = os.path.exists("/usr/bin/time") and os.uname().sysname == "Linux"

def run(cmd, env_extra=None):
    env = dict(os.environ)
    if env_extra: env.update(env_extra)
    if TIMEV: cmd = ["/usr/bin/time", "-v"] + cmd          # server protocol: peak RSS in stderr
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", p.stderr)
    rss_gb = int(m.group(1)) / 1048576.0 if m else 0.0
    return p.stdout, p.stderr, time.time() - t0, rss_gb

def dist_blocks(stdout):
    """split stdout into per-cell {core: count} dicts, in print order"""
    blocks, cur = [], None
    for ln in stdout.splitlines():
        if ln.startswith("[sct-peel] Max core:"):
            if cur is not None: blocks.append(cur)
            cur = {}
        m = re.match(r"core=(-?\d+) count=(\d+)", ln)
        if m and cur is not None:
            cur[int(m.group(1))] = int(m.group(2))
    if cur is not None: blocks.append(cur)
    return blocks

def main():
    if len(sys.argv) < 5:
        print("usage: nsi_sweep_gate.py <binary> <graph> <r> <s0> <smax>"); sys.exit(2)
    binp, graph, r, s0, smax = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    print(f"=== SWEEP {graph} r={r} s={s0}..{smax} ===", flush=True)
    so, se, tsweep, rss_sw = run([binp, graph, str(r), str(s0)], {"SCT_SWEEP": str(smax)})
    sweep_cells = dist_blocks(so)
    for ln in so.splitlines():
        if ln.startswith("[nsi") or "TIMING" in ln: print("  " + ln)
    if len(sweep_cells) != smax - s0 + 1:
        print(f"GATE FAIL: sweep produced {len(sweep_cells)} cells, expected {smax-s0+1}"); print(se[-3000:]); sys.exit(1)
    fails = 0
    tnative_sum = 0.0
    for i, cs in enumerate(range(s0, smax + 1)):
        no, ne, tn, rss_n = run([binp, graph, str(r), str(cs)])
        tnative_sum += tn
        nat = dist_blocks(no)
        natd = nat[-1] if nat else {}
        swd = dict(sweep_cells[i]); swd.pop(0, None)          # drop core-0 (absent-by-construction natively)
        natd.pop(0, None)
        if swd == natd:
            print(f"  cell s={cs}: GATE OK ({len(swd)} core levels, {sum(swd.values()):.0f} r-cliques)  native={tn:.2f}s rss={rss_n:.1f}GB")
        else:
            fails += 1
            allk = sorted(set(swd) | set(natd))
            bad = [(k, swd.get(k), natd.get(k)) for k in allk if swd.get(k) != natd.get(k)]
            print(f"  cell s={cs}: GATE FAIL  ({len(bad)} differing levels; first 8: {bad[:8]})")
    print(f"=== sweep-total={tsweep:.2f}s (rss={rss_sw:.1f}GB)  native-cells-total={tnative_sum:.2f}s  {'ALL GATES PASS' if fails==0 else str(fails)+' CELLS FAIL'} ===")
    sys.exit(1 if fails else 0)

if __name__ == "__main__":
    main()
