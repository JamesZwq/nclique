#!/usr/bin/env python3
# Multi-r NSI2 correctness gate: for each (r,s) in the plane, dump per-clique cores with the
# INDEPENDENT reference impl (build/bin/degeneracy_cliques, PIVOTER_RUN_REF), then query the
# serialized NSI2 index for every clique and require bit-exact agreement.
import subprocess, os, sys, random

REPO = "/Users/zhangwenqian/UNSW/pivoter"
REF  = f"{REPO}/build/bin/degeneracy_cliques"
NSIQ = f"{REPO}/region_native/nsi_query"
GR   = f"{REPO}/data/ca-GrQc.edges"
IDX  = "/tmp/grqc.nsi"
RMIN, RMAX, SMAX = 3, 5, 6
MAXQ = 100000

def load_ref(fn):
    d = {}
    for ln in open(fn):
        if ln.startswith('#'): continue
        p = ln.rstrip('\n').split('\t')
        if len(p) != 2: continue
        d[tuple(sorted(int(x) for x in p[0].split()))] = int(float(p[1]))
    return d

fails = tot = cells = 0
for r in range(RMIN, RMAX+1):
    for s in range(r+1, SMAX+1):
        ref = f"/tmp/nsi2_ref_{r}_{s}.core"
        env = dict(os.environ); env["PIVOTER_RUN_REF"]="1"; env["PIVOTER_DUMP_CORE"]=ref
        p = subprocess.run([REF, GR, str(r), str(s), "default"], capture_output=True, text=True, env=env)
        if not os.path.exists(ref):
            print(f"  r={r} s={s}: REF DUMP FAILED ({p.stderr.strip()[:80]})"); fails+=1; continue
        truth = load_ref(ref)
        keys = list(truth.keys()); random.seed(7); random.shuffle(keys); keys = keys[:MAXQ]
        if not keys:
            print(f"  r={r} s={s}: (no r-cliques)"); continue
        qf = f"/tmp/nsi2_q_{r}_{s}.txt"
        with open(qf,"w") as f:
            for k in keys: f.write(" ".join(map(str,k))+"\n")
        q = subprocess.run([NSIQ, IDX, "pointfile", str(r), str(s), qf], capture_output=True, text=True)
        got = [x for x in q.stdout.split()]
        if len(got) != len(keys):
            print(f"  r={r} s={s}: GATE FAIL (got {len(got)} answers != {len(keys)} queries; err={q.stderr.strip()[:80]})"); fails+=1; continue
        bad = sum(1 for k,g in zip(keys,got) if int(float(g)) != truth[k])
        tot += len(keys); cells += 1
        if bad:
            # show a couple mismatches
            ex = [(k,truth[k],got[i]) for i,(k,g) in enumerate(zip(keys,got)) if int(float(g))!=truth[k]][:3]
            print(f"  r={r} s={s}: GATE FAIL ({bad}/{len(keys)} mismatch) e.g. {ex}"); fails+=1
        else:
            print(f"  r={r} s={s}: BIT-EXACT ({len(keys)} cliques)")
print(f"=== {'ALL CELLS PASS' if fails==0 else str(fails)+' CELLS FAIL'} : {tot} point queries across {cells} cells ===")
sys.exit(1 if fails else 0)
