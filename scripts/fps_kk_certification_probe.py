#!/usr/bin/env python3
# FPS KK-sandwich certification test. For a fixed r and cell (r,s), use the EXACT per-r-clique cores
# at (r,s-1),(r,s),(r,s+1) (from PIVOTER_RUN_REF dumps, keyed by vertex set = s-invariant) to compute
# the Kruskal-Katona s-neighbor bounds and measure what fraction of r-cliques are CERTIFIED
# (f_lo == f_hi == kappa) WITHOUT peeling cell (r,s) -- i.e. bracketed exactly by the s-neighbors.
# g_a(m) = KK lower shadow (Lovasz form): C(x,a)=m (real x) -> C(x,a-1). Lemma 1 (SigmodPlus 128):
#   kappa_{r,s}(R)   >= g_{s+1-r}(kappa_{r,s+1}(R))         [lower from the s+1 neighbor]
#   kappa_{r,s-1}(R) >= g_{s-r}(kappa_{r,s}(R))  => upper: kappa_{r,s} <= g_{s-r}^{-1}(kappa_{r,s-1}(R))
import sys, math

def Creal(x, a):            # C(x,a) for real x>=0, integer a>=0
    if a == 0: return 1.0
    if x < a - 1e-9: return 0.0
    p = 1.0
    for i in range(a): p *= (x - i)
    return p / math.factorial(a)

def inv_C(m, a):            # smallest real x>=a-1 with C(x,a)=m  (a>=1, m>=0)
    if a == 0: return 0.0
    if m <= 0: return float(a - 1)
    lo, hi = float(a - 1), float(a)
    while Creal(hi, a) < m: hi *= 2
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if Creal(mid, a) < m: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def g(m, a):               # KK lower shadow from uniformity a to a-1
    if m <= 0: return 0.0
    if a <= 1: return m
    x = inv_C(m, a)
    return Creal(x, a - 1)

def load(fn):
    d = {}
    try: f = open(fn)
    except: return None
    for ln in f:
        if ln.startswith('#') or not ln.strip(): continue
        parts = ln.rstrip('\n').split('\t')
        if len(parts) != 2: continue
        vs = tuple(sorted(int(x) for x in parts[0].split()))
        d[vs] = float(parts[1])
    return d

def analyze(tag, r, s):
    km1 = load(f"/tmp/fps/{tag}_{r}_{s-1}.core")
    k0  = load(f"/tmp/fps/{tag}_{r}_{s}.core")
    kp1 = load(f"/tmp/fps/{tag}_{r}_{s+1}.core")
    if not (km1 and k0 and kp1):
        print(f"{tag} (r={r},s={s}): MISSING dumps"); return
    a_up, a_lo = s - r, s + 1 - r
    tot = cert = bracket_only_lo = bracket_only_hi = both = 0
    for R, kappa in k0.items():
        tot += 1
        # lower bound from s+1 neighbor (0 if R absent at s+1 -> no info)
        flo = g(kp1[R], a_lo) if R in kp1 else 0.0
        # upper bound from s-1 neighbor (+inf if R absent at s-1 -> no info)
        fhi = inv_and_up(km1[R], a_up) if R in km1 else float('inf')
        flo_i, fhi_i = round(flo), (round(fhi) if fhi != float('inf') else None)
        unsound = 0
        if fhi_i is not None and flo_i == fhi_i:
            cert += 1
            if kappa > 0.5: both += 1
            if abs(flo_i - round(kappa)) > 0: bracket_only_lo += 1   # SOUNDNESS: certified must equal true kappa
    nz = sum(1 for v in k0.values() if v>0.5)
    print(f"{tag} cell (r={r},s={s}): r-cliques={tot} (nz-core={nz})  KK CERTIFIED={100.0*cert/tot:.1f}% ({cert}/{tot}) | nz-core certified={100.0*both/max(nz,1):.1f}% | UNSOUND(cert!=true kappa)={bracket_only_lo}")

def inv_and_up(m, a):      # upper bound: g_a^{-1}(m) = C(x,a) where C(x,a-1)=m
    if m <= 0: return 0.0
    if a <= 1: return m
    x = inv_C(m, a - 1)
    return Creal(x, a)

if __name__ == "__main__":
    import sys
    tags = sys.argv[1].split(",") if len(sys.argv)>1 else ("grqc","hepph","epin")
    for tag in tags:
        analyze(tag, 3, 5)
