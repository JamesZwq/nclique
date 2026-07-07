#!/usr/bin/env python3
# §136 NSI index gate: build the index once (SCT_INDEX_OUT + SCT_SWEEP), REF-dump the exact
# per-r-clique cores for every cell, then query the index for EVERY dumped r-clique at EVERY
# cell and require exact agreement. Also runs the latency bench and prints the stats line.
#   python3 nsi_index_gate.py <sweep_binary> <main_binary> <graph> <r> <s0> <smax> [maxq]
import subprocess, sys, os, random, time

sweepbin, mainbin, graph, r, s0, smax = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
maxq = int(sys.argv[7]) if len(sys.argv) > 7 else 2000000
NSIQ = os.path.join(os.path.dirname(sweepbin) or ".", "nsi_query")
idx = f"/tmp/nsi_{os.path.basename(graph)}_{r}_{s0}_{smax}.nsi"

print(f"=== build index: {graph} r={r} s={s0}..{smax} ===", flush=True)
env = dict(os.environ); env["SCT_SWEEP"] = str(smax); env["SCT_INDEX_OUT"] = idx
p = subprocess.run([sweepbin, graph, str(r), str(s0)], capture_output=True, text=True, env=env)
for ln in p.stdout.splitlines():
    if "[nsi-index]" in ln or "[nsi-sweep]" in ln: print("  " + ln)
if p.returncode != 0: print("BUILD FAILED"); sys.exit(1)

def load(fn):
    d = {}
    for ln in open(fn):
        if ln.startswith('#'): continue
        parts = ln.rstrip('\n').split('\t')
        if len(parts) != 2: continue
        d[tuple(sorted(int(x) for x in parts[0].split()))] = int(float(parts[1]))
    return d

fails = tot = 0
qfile = "/tmp/nsi_queries.txt"
for s in range(s0, smax + 1):
    ref = f"/tmp/nsi_ref_{s}.core"
    env2 = dict(os.environ); env2["PIVOTER_RUN_REF"] = "1"; env2["PIVOTER_DUMP_CORE"] = ref
    subprocess.run([mainbin, graph, str(r), str(s), "default"], capture_output=True, text=True, env=env2)  # "default" sort keeps ORIGINAL vertex ids (degen sort relabels)
    truth = load(ref)
    keys = list(truth.keys())
    random.seed(42); random.shuffle(keys)
    keys = keys[:maxq]
    with open(qfile, "w") as f:
        for k in keys: f.write(" ".join(map(str, k)) + "\n")
    # batch-query via the bench path is aggregate; for the GATE use a dedicated compare mode:
    # feed queries and expected cores to nsi_query via a python-side one-shot subprocess per CELL
    # would be slow; instead reuse 'bench' for latency and do correctness via 'spectrumfile'...
    # v1: correctness via per-cell batch: run nsi_query in 'pointfile' mode
    p = subprocess.run([NSIQ, idx, "pointfile", str(s), qfile], capture_output=True, text=True)
    got = [float(x) for x in p.stdout.split()]
    if len(got) != len(keys):
        print(f"  cell s={s}: GATE FAIL (answer count {len(got)} != {len(keys)})"); fails += 1; continue
    bad = sum(1 for k, g in zip(keys, got) if int(g) != truth[k])
    tot += len(keys)
    if bad: print(f"  cell s={s}: GATE FAIL ({bad}/{len(keys)} mismatches)"); fails += 1
    else:   print(f"  cell s={s}: GATE OK ({len(keys)} r-cliques exact)", flush=True)

print("=== stats ===")
subprocess.run([NSIQ, idx, "stats"])
print("=== latency (queries from the deepest cell's r-cliques) ===")
subprocess.run([NSIQ, idx, "bench", qfile])
print(f"=== {'ALL GATES PASS' if fails == 0 else str(fails) + ' CELLS FAIL'} ({tot} point queries checked) ===")
sys.exit(1 if fails else 0)
