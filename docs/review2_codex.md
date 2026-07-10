# Hostile Review #2 (codex gpt-5.6-sol, max) — 2026-07-10

**VERDICT: Reject, 2/5.** Credible algorithmic core (profile quotient, weighted-pattern peeling, shell-order
replay), but the central cross-r queryable-plane claim is unsupported: plane table all TBD, NSI2 unimplemented,
fresh results mixed runtime + worse memory, one theorem false at zero.

## Reject-triggering (ranked)
1. **"competitive" plane claim unsupported.** Table plane all TBD; fresh range 0.38x–3.06x (dense wins, sparse
   loses) + worse memory (447GB vs 401GB HepPh; 5.5 vs 2.6GB dblp). Not evidence for unqualified "competitive."
   FIX: populate table w/ exact ratios+absolute times+RSS, 3 trials, define FixedRows-mem aggregation, phase
   attribution, characterize dense-wins/sparse-loses, replace "competitive" w/ measured range.
2. **Queryable plane index (NSI2) doesn't exist.** Query numbers are warm fixed-r-column lookup-kernel only
   (200k pre-validated cliques, excludes validation/alloc/deser/load/IO). Intro/abstract claim "any (r,s) query
   in ns." FIX: implement NSI2 (serialize shared header + all r-columns, random-r correctness, total bytes,
   load time, warm/cold, validation-inclusive median+p95, 3+ trials) OR reframe as "fixed-row decomposition +
   proposed plane layout" and delete plane-query claims.
3. **Theorem FALSE at zero. [VERIFIED REAL]** g_a(0)=C(a-1,a-1)=1, so Shadow Bound asserts κ_{r,s-1}≥1 when
   κ_{r,s}=0 (K_r, s≥r+2: both zero → 0≥1 false). Also Shadow/Diagonal proofs omit global feasibility
   construction; missing from formal_theory.tex. FIX: define g_a(0)=0 (or restrict positive cores); prove Shadow
   via shadows of ALL witnesses; Diagonal via B=∪_R K_{r-1}(R); move proofs to supplement. Chain salvageable.
4. **Novelty not established.** Construction = "weighted-class version of canonical Pivoter recurrence";
   KClist++ undiscussed; no quotient-disabled ablation. FIX: theorem-by-theorem prior-art table (nucleus,
   Pivoter/SCT, KClist++, core/truss indexes) + quotient-off measurements (classes, raw/compressed leaves,
   cliques/patterns, build time, RSS, residue).

## Major
1. FixedRows comparison empty → the two-baseline framing looks defensive; 1600x validates fixed-row engine not
   cross-r structure; add phase + shared-preprocessing controls.
2. CND 96-core vs our 1-thread: doesn't weaken elapsed-time (conservative) but invalidates clean mem/feasibility
   attribution (96-thread buffers may explain part of 140x + 300GB fails). Report CND @1/mid/96 threads + CPU-sec.
3. Submission visibly unfinished (compiler/cmds/seeds TODO; -march=native fresh vs -O3 in paper; no trials/variance;
   plane gate TODO). [NOTE: these are LaTeX COMMENTS invisible in PDF; the empty TABLE is the real visible flaw.]
4. Plane evidence absent from paper; fresh covers 3 graphs r=3..5 only; no ablation isolating class compression.
5. Cost analysis tautological (T_shared+ΣT_col); no bound/count for regions/profiles/leaves/patterns/antichain.
   447GB shows it's not cosmetic. Add parameterized complexity + per-phase time/RSS + failure thresholds.
6. Size bound Eq undercounts: O(|L|),O(|P_r|) ignore per-leaf class intervals + per-pattern variable composition.
   Use total leaf-class incidences + pattern nonzeros.
7. Index-size inconsistent (abstract 0.001–2 vs Exp 0.001–2.2); denominator is fixed-r column not plane.
8. Hierarchy stores scalar cores not connected-nucleus IDs; comparison w/ component-query core/truss indexes
   incomplete. Implement connected-nucleus retrieval OR restrict claims to core-value lookup.

## Minor
Cross-cell nesting g noninteger (use ceil, k≥1); Diagonal restrict r≥2; Construction op contracts missing;
budget 300GB vs CND "requires 310GB" (mark censored); CaseStudy II 54% vs 98% mismatched cohorts; CaseStudy I
Jaccard needs tie-break; "13 graphs + ca-GrQc = 14" state plainly.

## Strongest / protect
The exact quotient-peeling argument (NOT the inherited pivot recurrence): weighted-pattern peeling equivalence
(uniqueness of largest k-feasible set → superlevel sets are unions of pattern orbits), class-SCT exact-cover,
shell-order replay invariants. Preserve these + honest limitations + two-baseline separation + fixed-row ablation.

## What flips to accept
1. Implement+validate NSI2 (real serialized queryable multi-r index, full size + end-to-end query).
2. Replace TBD plane table w/ multi-trial results, report 0.38x–3.06x range + worse-mem + phase attribution + gates.
3. Repair false zero-case shadow, complete Shadow/Diagonal feasible-set proofs, fix hierarchy threshold + size bound.
4. Establish novelty vs Pivoter/SCT/KClist++/nucleus/FP-indexes + quotient-off ablation + compression stats.
5. Normalize CND across thread counts, report CPU work + memory-feasible config + exact cmds/versions/raw/variance.
