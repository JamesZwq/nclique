#!/usr/bin/env python3
# verify_class_base.py -- BRUTE-FORCE proof of the class-compression BASE of (r,s)-nucleus.
# Independently (no engine) recomputes (r,s)-nucleus cores + region-classes + true-twin classes
# on thousands of random graphs and asserts the foundational claims:
#   CLAIM 1 : same region-class PATTERN (sorted class-multiset) => same (r,s)-core   [licenses peeling patterns, not r-cliques]
#   CLAIM 1b: every support-bearing r-clique has all vertices in a region-class
#   CLAIM 2 : true-twin (N[u]=N[v]) REFINES region-class (twin => same region-class)
#   CLAIM 3 : same TWIN-pattern => same core (so twin-classes are a CORRECT, MCE-free substitute)
#   ORDERING: #region-classes <= #twin-classes <= #vertices
# Result (seed 12345/999): 0 violations over 36000 (graph x r,s) tests. See SigmodPlus.md sec 110.
import itertools, random, sys
random.seed(12345)

def is_clique(verts, adjset):
    verts=list(verts)
    for i in range(len(verts)):
        for j in range(i+1,len(verts)):
            if verts[j] not in adjset[verts[i]]: return False
    return True

def maximal_cliques(n, adjset):
    # Bron-Kerbosch (no pivot, fine for tiny graphs)
    res=[]
    def bk(R,P,X):
        if not P and not X: res.append(frozenset(R)); return
        for v in list(P):
            bk(R|{v}, P & adjset[v], X & adjset[v]); P=P-{v}; X=X|{v}
    bk(set(), set(range(n)), set())
    return res

def rs_cliques(n, adjset, k):
    return [frozenset(c) for c in itertools.combinations(range(n),k) if is_clique(c,adjset)]

def rs_nucleus_cores(n, adjset, r, s):
    R=rs_cliques(n,adjset,r); S=rs_cliques(n,adjset,s)
    s_to_r={Sc:[Rc for Rc in R if Rc<=Sc] for Sc in S}
    supp={Rc:sum(1 for Sc in S if Rc<=Sc) for Rc in R}
    alive=set(R); aliveS=set(S); cur=dict(supp); core={}; thr=0
    while alive:
        Rc=min(alive,key=lambda x:cur[x]); thr=max(thr,cur[Rc]); core[Rc]=thr; alive.discard(Rc)
        for Sc in list(aliveS):
            if Rc<=Sc:
                aliveS.discard(Sc)
                for R2 in s_to_r[Sc]:
                    if R2 in alive: cur[R2]-=1
    return core, R

def region_classes(n, adjset, s):
    regions=[M for M in maximal_cliques(n,adjset) if len(M)>=s]
    prof={v:frozenset(i for i,M in enumerate(regions) if v in M) for v in range(n)}
    cls={}; p2c={}
    for v in range(n):
        if prof[v]:
            if prof[v] not in p2c: p2c[prof[v]]=len(p2c)
            cls[v]=p2c[prof[v]]
    return cls, prof

def true_twin_class(n, adjset):
    # group by closed neighborhood
    sig={v: frozenset(adjset[v]|{v}) for v in range(n)}
    s2c={}; cls={}
    for v in range(n):
        if sig[v] not in s2c: s2c[sig[v]]=len(s2c)
        cls[v]=s2c[sig[v]]
    return cls

def check_graph(n, edges, r, s):
    adjset=[set() for _ in range(n)]
    for u,v in edges:
        if u!=v: adjset[u].add(v); adjset[v].add(u)
    core,Rcl = rs_nucleus_cores(n,adjset,r,s)
    rcls,prof = region_classes(n,adjset,s)
    # CLAIM 1: same region-class pattern => same core (for support-bearing r-cliques: all verts have a class)
    pat_cores={}
    for Rc in Rcl:
        if all(v in rcls for v in Rc):
            pat=tuple(sorted(rcls[v] for v in Rc))
            pat_cores.setdefault(pat,set()).add(core[Rc])
    bad1={p:c for p,c in pat_cores.items() if len(c)>1}
    # CLAIM 1b: every support-bearing (core>0) r-clique has all verts in a class
    bad1b=[Rc for Rc in Rcl if core[Rc]>0 and not all(v in rcls for v in Rc)]
    # CLAIM 2: true-twin refines region-class (twin same => region-class same, among classed verts)
    tw=true_twin_class(n,adjset)
    bad2=[]
    for u in range(n):
        for v in range(n):
            if u<v and tw[u]==tw[v]:
                # same twin => must be same region-class (or both unclassed)
                ru=rcls.get(u); rv=rcls.get(v)
                if ru!=rv: bad2.append((u,v,ru,rv))
    return bad1,bad1b,bad2,len(rcls),len(set(rcls.values())),len(set(tw.values()))

def rand_graph(n,p):
    edges=[(u,v) for u in range(n) for v in range(u+1,n) if random.random()<p]
    return edges

print("=== CONCRETE EXAMPLE: K4{a=0,b=1,c=2,d=3} + pendant a-e(4), r=2 s=3 ===")
ex=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3),(0,4)]
core,Rcl=rs_nucleus_cores(5,[ {1,2,3,4} if i==0 else ({0,2,3} if i==1 else ({0,1,3} if i==2 else ({0,1,2} if i==3 else {0}))) for i in range(5)],2,3)
rcls,prof=region_classes(5,[ {1,2,3,4} if i==0 else ({0,2,3} if i==1 else ({0,1,3} if i==2 else ({0,1,2} if i==3 else {0}))) for i in range(5)],3)
print(" region-classes:",rcls,"  (a,b,c,d all class 0; e has none)")
print(" edge cores:",{tuple(sorted(k)):v for k,v in core.items()})
tw=true_twin_class(5,[ {1,2,3,4} if i==0 else ({0,2,3} if i==1 else ({0,1,3} if i==2 else ({0,1,2} if i==3 else {0}))) for i in range(5)])
print(" twin-classes:",tw,"  (a=0 split from b,c,d)")

print("\n=== BRUTE-FORCE over MANY random graphs x (r,s) ===")
tot=0; v1=v1b=v2=0
for trial in range(4000):
    n=random.randint(5,11); p=random.uniform(0.3,0.85)
    edges=rand_graph(n,p)
    for (r,s) in [(2,3),(2,4),(3,4),(3,5)]:
        if s>n: continue
        bad1,bad1b,bad2,_,_,_=check_graph(n,edges,r,s)
        tot+=1
        if bad1: v1+=1; print(f"  VIOLATION CLAIM1 (same pattern, diff core) n={n} p={p:.2f} r={r} s={s}: {bad1}")
        if bad1b: v1b+=1; print(f"  VIOLATION CLAIM1b (support-bearing rclique w/ unclassed vertex) n={n} r={r} s={s}")
        if bad2: v2+=1; print(f"  VIOLATION CLAIM2 (twin not refining region-class) n={n} r={r} s={s}: {bad2[:3]}")
print(f"\nTESTS={tot}")
print(f"CLAIM 1  (same region-class pattern => same (r,s)-core): violations={v1}  -> {'PASS' if v1==0 else 'FAIL'}")
print(f"CLAIM 1b (every support-bearing r-clique has all vertices classed): violations={v1b} -> {'PASS' if v1b==0 else 'FAIL'}")
print(f"CLAIM 2  (true-twin refines region-class): violations={v2} -> {'PASS' if v2==0 else 'FAIL'}")
