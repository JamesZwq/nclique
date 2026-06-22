#!/usr/bin/env python3
# Validate direction (2): can a pivot/hold bundle's peel output be compressed?
#  A: class compression  = #patterns / #r-cliques (sum mult).
#  B: per-bundle OUTPUT compression = Σ#distinct-core / Σ#patterns (<<1 = headroom).
#  C: within-bundle core monotone in support? (bands findable by support order)
#  D: within-bundle core a FUNCTION of support? (#distinct(support,core)==#distinct support)
#     -> the killer: if yes, the whole bundle's cores come from the support formula.
# Bundle = patterns sharing reg0 (min host region). Reported for ALL and for hostSz==1 (clean).
import sys
from collections import defaultdict

path = sys.argv[1]
npat = 0; rcliques = 0; multi_host = 0
bundles_all = defaultdict(list)   # reg0 -> [(sup,core)]
bundles_sh  = defaultdict(list)   # reg0 -> [(sup,core)] for hostSz==1

with open(path) as f:
    f.readline()
    for line in f:
        p = line.split('\t')
        if len(p) < 6: continue
        reg0=int(p[0]); hostSz=int(p[1]); sup=int(p[2]); core=int(float(p[3])); mult=int(p[4])
        npat += 1; rcliques += mult
        bundles_all[reg0].append((sup,core))
        if hostSz==1: bundles_sh[reg0].append((sup,core))
        else: multi_host += 1

print(f"=== {path} ===")
print(f"#patterns          = {npat:,}")
print(f"#r-cliques (Σmult) = {rcliques:,}")
print(f"class compression  = {npat/max(1,rcliques):.4f}  (patterns/r-cliques; <1=compression)")
print(f"multi-host         = {multi_host:,} ({100*multi_host/max(1,npat):.1f}%)")

def analyze(bundles, label):
    tot_pat=0; tot_core=0
    big=[]; mono_ok=0; mono_bad=0; func_ok=0; func_bad=0; multi_core_bundles=0
    for reg0,lst in bundles.items():
        n=len(lst); tot_pat+=n
        cores=set(c for _,c in lst); nd=len(cores); tot_core+=nd
        if n>=20: big.append((reg0,n,nd))
        if nd>=2:
            multi_core_bundles+=1
            # D: core a function of support?
            sc=set(lst); ns=len(set(s for s,_ in lst))
            if len(sc)==ns: func_ok+=1
            else: func_bad+=1
            # C: core monotone non-decreasing in support? (O(n log n))
            minc=defaultdict(lambda:10**9); maxc=defaultdict(lambda:-1)
            for s,c in lst:
                if c<minc[s]: minc[s]=c
                if c>maxc[s]: maxc[s]=c
            run=-1; ok=True
            for s in sorted(minc):
                if minc[s]<run: ok=False; break
                run=max(run,maxc[s])
            if ok: mono_ok+=1
            else: mono_bad+=1
    print(f"\n-- {label}: bundles={len(bundles):,} Σpat={tot_pat:,} Σdistinct-core={tot_core:,}")
    if tot_pat: print(f"   OUTPUT compression = {tot_core/tot_pat:.4f} (Σdistinct-core/Σpat; <<1=direction-2 headroom)")
    if big:
        bp=sum(n for _,n,_ in big); bd=sum(d for _,_,d in big)
        print(f"   big bundles(>=20)  = {len(big):,} Σpat={bp:,} Σdistinct={bd:,} ratio={bd/bp:.4f}")
        big.sort(key=lambda x:-x[1])
        for reg0,n,d in big[:6]: print(f"      reg {reg0}: {n} pat -> {d} cores ({d/n:.3f})")
    if multi_core_bundles:
        print(f"   multi-core bundles = {multi_core_bundles:,}")
        print(f"   C monotone(core↑sup): ok={mono_ok:,} bad={mono_bad:,}  frac={mono_ok/multi_core_bundles:.3f}")
        print(f"   D core=f(support):    ok={func_ok:,} bad={func_bad:,}  frac={func_ok/multi_core_bundles:.3f}")

analyze(bundles_sh,  "Metric B/C/D: single-host bundles (hostSz==1, clean)")
analyze(bundles_all, "Metric B/C/D: ALL patterns by min-region (coarse)")
