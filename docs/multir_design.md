# Multi-r Whole-Plane Index -- codex (gpt-5.6-sol) design (2026-07-10)

## Verdict
- The class-SCT LEAF representation is r-INDEPENDENT: count_{L,s}(P) = sum over y (sum y_c = s,
  max(l_c,b_c) <= y_c <= u_c) of prod_c C(n_{L,c}-b_c, y_c-b_c). Pattern supplies b (|b|=r), leaf
  supplies n,l,u, cell supplies T=s. This is the DP in CCPathCore.h:125 (count_with_extra_lower_impl),
  summed in sctSupport (region_native_sct_peel.cpp:1336). NO r stored in a leaf.
- The builder scalableBuildClassSCT(C,weights,adj,k,kOver) takes NO r -- only clique-size range
  [k,kOver] (ClassSCTScalable.h:314). One tree serves all s in [k,kOver] (already exploited across s:
  baseLeaves built once @964, retargeted per cell @1450).
- BUT the current executable is r-dependent: it does r-mergeable split (@425-576) and REMOVES those
  regions (@517) BEFORE building classes/quotient/SCT (@578-617, @945-973). So the built tree differs
  per r. FIX: build shared class map + tree from the UNSPLIT region set FIRST, then per-r views.
- The paper's `\binom{w}{s-r}` is only the SPECIAL CASE (single unrestricted pool); the general leaf
  count is the product-of-binomials above. State the general form first (Construction.tex:72-75).
- BuildCPI's `order r` input is spurious (never used) -> should take [S_min,S_max].

## Genuinely r-dependent state (stays per-r)
r-mergeable classification (|M cap M'| < r); pattern set P_r (compositions sum to r); multiplicity
prod C(w_c,b_c); pattern-to-leaf maps; c(P); peel/forbidden state (witnessTail=s-r). Optional amortize:
mu(M)=max_{M'} |M cap M'|; M is r-mergeable iff mu(M) < r (pairwise fallback already computes it @503
but discards it).

## Diagonal theorem = UPPER bound only
U_rho(Q) = min_{c:q_c>0} kappa_{rho-1,rho}(Q - e_c) - 1  (TransferTheory.tex:188).
L_rho(Q) = c(Q) - rho (clique floor). If U=L: boundary exact = L (certify). If U>L: must peel.
U<L: invariant failure (bug). Do NOT schedule at U (upper bound != death level) -- not justified.

## Multi-r whole-plane pseudocode (BuildPlaneNSI)
Smin = rmin+1; Regions = maximal cliques >= Smin; (classOf,weights,profiles,Q)=BuildClasses(Regions);
SharedTree = BuildClassSCT(Q,weights,Smin,Smax)  // immutable, no r-patterns/peel state.
for rho = rmin..rmax:
  boundary = rho+1; Mergeable_rho = classify(Regions,rho); Patterns_rho = enumerate(Regions,weights,rho)
  partition -> Direct_rho (unique mergeable host), Active_rho; c[P]=max{|M|:M hosts P}; maps_rho over SharedTree
  if rho=rmin: Kbnd = FullPeel(boundary)
  else:
    for Q in Active_rho: U=min_c (LookupExact(PrevDiagonal, rho-1, Q-e_c) - 1); L=c[Q]-rho; assert U>=L
      if U=L: Kbnd[Q]=L; Certified.add(Q,death=L)  else Residue.add(Q)
    Kbnd[Residue] = ReplayPeel((rho,boundary), Residue, scheduledExactDeaths=Certified)
  direct mergeable: Kbnd[P]=|M|-rho; closedFrom[P]= boundary if Kbnd=c[P]-rho else +inf
  PrevDiagonal = exact view of Kbnd (incl mergeable closed forms); Prev=Kbnd
  for s=boundary+1..Smax:
    for P in Active_rho: if Prev[P]=C(c[P]-rho,s-1-rho): K[P]=C(c[P]-rho,s-rho); Certified.add; closedFrom
       else Residue.add
    K[Residue]=ReplayPeel((rho,s),Residue,Certified); store only residue K; Prev=K
  emit Column_rho(pattern table, boundary, closedFrom, residue dicts)

## Storage NSI2 + query
Shared: classOf, weights, class->region profile CSR, region sizes, optional mu(M), ONE shared tree
(counted once). Per-r column: boundary b_r=r+1, pattern dict comp->id, per active P {comp,mult,c(P),
boundary core, closedFrom}, residue dict D_{r,s}. Point query O(r+profile-intersect); with closedFrom
-> arithmetic or one residue lookup. Active-pattern and mergeable cases disjoint.

## Cost (HONEST -- paper-safe claim)
T = MCE + classes + shared tree + sum_r ( merge_r + patterns_r + maps_r + boundary_r +
    sum_{s>b_r}(O(|P_r|)+residue_peel) ). NOT unconditionally "one cell per column": only when residues
empty. Worst case (no cert) = full peel every cell; |P_r| can be combinatorial. CENTRAL CLAIM:
"One universal CPI is a sufficient counting representation for the entire bounded (r,s) plane." The
"plane is cheap" part must be shown EXPERIMENTALLY (diag+chain cert rates, residue sizes, time, memory).

## 12-step region_native_sct_peel.cpp refactor
1. Parse rMin,rMax,Smax,Smin=rMin+1 (keep fixed-r as compat wrapper). @360
2. MCE once at Smin (change @412-419 threshold from boundary s to Smin); regions immutable.
3. Move class/profile construction (@578-617) BEFORE r-mergeability; r-merge (@425-576) must NOT
   remove regions from the global vector.
4. Build quotient+tree ONCE before the r-loop (@945-973): scalableBuildClassSCT(nC,weights,qadj,Smin,Smax)
   -> immutable sharedLeaves.
5. r-mergeability = per-column mask: MergeView classifyMergeable(SharedRegions&, r) returning flags +
   direct metadata, no mutation. Keep SCT_NO_RMERGE ablation. Optional maxOverlap[M] path.
6. buildPatterns(r, mergeView) wrapping @778-914: hosts, multiplicity, c(P)=max_{M in host}|M|.
7. buildPatternLeafMaps(sharedLeaves,patterns,r,r+1,Smax) @1005-1421; per-r runtime leaf VIEW, strip
   metadata only from the view (do NOT strip sharedLeaves @1268-1272).
8. Nest current s-loop (@1450) inside an outer r-loop; per (r,s) reset pattern/peel state (deadY fresh);
   keep nonmutating a_Y path (@1763).
9. Generalize fpsCell (@1460-1494) to a seeded replay cell: boundary branch for r>rmin computes U,L from
   PrevDiagonal, certifies U=L, else residue; reuse skip-init @1563-1568, scheduled deaths @1480-1489,
   queue omit @2135-2149, leaf skip @2544-2548.
10. Keep PrevDiagonal composition->core view alive between columns (active + mergeable closed forms).
11. Record closedFrom during the sweep (first cell where exact core = clique floor).
12. NSI2 serialization (replace NSI1 @3458-3484): shared header + shared classes/profiles/regions +
    optional shared tree + #r-columns + directory of per-r offsets + per-r blocks.

## CORRECTNESS GATE (definition of done)
Plane mode with rMin=rMax=R must be BIT-EXACT vs the current fixed-r engine at r=R (same core-value
distribution per cell) on >=3 gate graphs (e.g. ca-GrQc, ca-HepPh, com-dblp). And a plane run over
r=3..5 must, per r-column, match the standalone fixed-r run at that r.
