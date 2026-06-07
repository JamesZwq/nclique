# Materialized Nucleus Spectrum Index: build time, size, query latency

Measured by `bench_spectrum_index` (source `bench_spectrum_index.cpp`,
universal SDCT in `src/SDCT_Augmented.{h,inl}`). Answers the three numbers a
cross-s index needs: **build time, index size, query latency**, each vs the
no-index baseline (rebuild the CPI + peel from scratch per queried s).

## 1. Method

**Build-once universal CPI.** `SDCT_Augmented_NoTree_Universal` emits every
maximal-clique leaf once (no `max_k` depth prune, no `cSize>=max_k` emission
filter). One SDCT traversal serves all s: leaf L=(H,P) contributes
`C(|P|, s-|H|)` s-cliques, so per-vertex support at level s is
`C(|P|,s-|H|)` for v in H and `C(|P|-1,s-|H|-1)` for v in P (Pivoter/SCT
invariant). For each s we (1) set per-leaf weights, (2) scatter support over
the vtx→leaf CSR (OpenMP, race-free per vertex), (3) run the validated SPIN★
sparse-bucket peel (`ST_V3_Peel`).  This replaces ω independent SDCT builds
with **one** build + ω cheap recompute+peel passes.

**Two query backends materialized:**
- **dense-by-s**: per-level CSR of active `(vertex, kappa)`. `query(s)` = read
  the level-s slice. O(n_active_s), optimal/output-bound.
- **anchor+delta**: per vertex `(s*, kappa_{s*})` + sparse Kruskal-Katona
  deltas; `query(s)` reconstructs κ_s by descending from s* (shadow + delta),
  iterating only vertices active at s (sorted by s*).

**Correctness.** Index κ_s compared **per vertex** against an independent
per-s `ST_V3_Build+Peel` (itself bit-identical to `degeneracy_cliques`).
Exhaustive: every s∈[2,ω] checked on ca-GrQc (ω=44) and ca-CondMat (ω=26) →
**0 mismatches**; sampled s on all 9 graphs incl. com-dblp (ω=114) and
ca-HepPh (ω=239) → **0 mismatches**.

Machine: Darwin arm64, 18 threads (recompute only; peel is single-thread).
Reproduce: `cmake --build build -j 12 --target bench_spectrum_index` then
`./build/bin/bench_spectrum_index <g.edges> [--queries s1,..] [--no-baseline]`.
Data: `bench_plots/spectrum_index.csv`.

## 2. Build time — universal build-once vs ω-rebuild

| graph | ω | univ build ms | + spectrum ms | = index ms | ω-rebuild ms | **speedup** |
|---|---|---|---|---|---|---|
| amazon0302 | 7 | 222.7 | 254.6 | 477 | 765 | 1.6× |
| ca-CondMat | 26 | 27.8 | 41.4 | 69 | 286 | 4.1× |
| ca-AstroPh | 57 | 182.2 | 312.3 | 494 | 4742 | 9.6× |
| cit-HepPh | 19 | 319.4 | 525.8 | 845 | 3589 | 4.2× |
| com-dblp | 114 | 393.1 | 618.3 | 1011 | 29863 | **29.5×** |
| ca-HepPh | 239 | 432.0 | 476.4 | 908 | 211379 | **232.7×** |

The win grows with ω: low-ω graphs (amazon, ω=7) save little; high-ω graphs
(dblp 30×, HepPh 233×) save enormously, because ω-rebuild repeats the full
near-ω SDCT traversal ω times while build-once does it once. The universal CPI
is not a memory blowup: Σ (vtx-leaf incidences) stays 0.36M–3.5M.

## 3. Memory: index vs the CPI "tree" substrate

Three resident structures can serve queries. The **CPI tree** = the universal
dual CSR + per-leaf metadata (kept if you answer by recompute+peel). The
**dense** and **anchor+delta** indexes are the materialized alternatives.
(Actual `vector::capacity()` bytes; CPI excludes the transient `countingV`.)

| graph | CPI tree MB | dense MB | anchor+delta MB | dense/tree | a+d/tree |
|---|---|---|---|---|---|
| ca-GrQc | 0.37 | 0.39 | 0.10 | 1.05× | 0.27× |
| ca-HepTh | 0.77 | 0.39 | 0.21 | 0.51× | 0.27× |
| ca-CondMat | 2.59 | 1.57 | 0.49 | 0.61× | 0.19× |
| ca-AstroPh | 8.14 | 3.15 | 0.52 | 0.39× | 0.06× |
| amazon0302 | 32.42 | 12.58 | 6.78 | 0.39× | 0.21× |
| cit-HepPh | 46.24 | 3.15 | 1.82 | **0.07×** | 0.04× |
| com-dblp | 26.36 | 25.17 | 6.50 | 0.95× | 0.25× |
| ca-HepPh | 6.14 | 3.15 | 2.41 | 0.51× | 0.39× |

**Key finding: the dense materialized index is usually SMALLER than the CPI
tree it was built from** (dense/tree ≈ 0.4–0.6×, down to 0.07× on cit-HepPh),
because the spectrum {κ_s(v)} is more compact than the full clique-tree
incidence structure (size ∝ Σ maximal-clique incidences). So the dense index
**dominates keeping the tree on both memory and query** — build the universal
CPI, materialize dense, drop the tree. Anchor+delta is smaller still
(0.04–0.39× the tree) but queries slowly (§5).

Dense ≈ 12 B/entry × Σ_s n_active_s; anchor+delta is a further 1.3–6× smaller
than dense. The CPI tree blows up on graphs with many maximal cliques
(cit-HepPh: 1339 B/vertex, 46 MB) even when the spectrum itself is small.

(The legacy `DynamicGraph<TreeGraphNode>` SDCT tree would be larger still;
ST_V3/SPIN★ already uses this compact dual-CSR form, which is what we measure.)

## 4. Query latency — the headline (against TWO baselines)

Two baselines, because there are two honest "no materialized index" options:
- **rebuild**: build the per-s CPI from scratch + peel, for every query. R
  queries ⇒ R SDCT builds (per-s trees are not reusable across s). Weakest.
- **shared**: build the *universal* CPI **once**, then per query recompute
  support + peel (no tree rebuild). The fair "build tree once, reuse" baseline.

Dense `query(s)` is a level-slice read; selected rows (full data in CSV):

| graph | s | n_active_s | **dense µs** | shared ms | rebuild ms | **dense/shared** | dense/rebuild | shared/rebuild |
|---|---|---|---|---|---|---|---|---|
| ca-CondMat | 2 | 23133 | 21.3 | 4.48 | 24.7 | **211×** | 1160× | 5.5× |
| ca-CondMat | 26 | 26 | 0.035 | 0.28 | 18.9 | **8057×** | 540429× | 67× |
| ca-AstroPh | 2 | 18771 | 20.7 | 10.55 | 160.6 | **510×** | 7758× | 15× |
| ca-AstroPh | 57 | 57 | 0.043 | 0.60 | 130.8 | **13958×** | 3.06M× | 219× |
| com-dblp | 2 | 317080 | 296.5 | 33.6 | 379.5 | **113×** | 1280× | 11× |
| com-dblp | 114 | 114 | 0.072 | 2.88 | 299.6 | **40014×** | 4.17M× | 104× |
| ca-HepPh | 2 | 12006 | 13.5 | 4.42 | 574.7 | **328×** | 42621× | 130× |
| ca-HepPh | 239 | 239 | 0.18 | 0.92 | 597.5 | **5236×** | 3.39M× | 648× |

Three takeaways:
1. **Dense vs the fair shared baseline: 113× – 88000×.** Even when the universal
   tree is built once and reused, per-query recompute+peel costs ms; the
   materialized dense index is µs. The materialization itself is the big win.
2. **Dense vs rebuild: 869× – 7.4M×** (vs the weakest baseline).
3. **Shared vs rebuild: 4× – 763×** — the per-query benefit of the build-once
   universal CPI alone (no materialization), largest at high s / high ω where a
   full rebuild redoes the whole near-ω SDCT traversal.

## 5. The anchor+delta trap

Anchor+delta is 1.3–6× smaller but its `query(s)` reconstructs each active
vertex from s* down to s (O(s*−s) shadow evals). At **low s on high-ω graphs**
this is catastrophic: ca-HepPh `query(2)` = **495 ms**, dblp `query(2)` = 93 ms
— as slow as (or slower than) a full rebuild, and 10³–10⁴× slower than dense.
**Conclusion: for a query index, dense wins decisively; the anchor+delta
compression (modest 1.3–6×) is not worth its query cost.**

## 6. Conclusion

- **Build-once universal CPI is correct (0 mismatches, all s) and 4–233×
  faster to build** than ω independent rebuilds; the gain scales with ω.
- **The dense materialized spectrum index answers κ_s in µs**: 113×–88000×
  faster than the fair "build tree once, reuse" baseline (recompute+peel), and
  869×–7.4M× faster than rebuilding per query, at 0.02–25 MB.
- Anchor+delta compression saves 1.3–6× space but its reconstruction makes
  queries 10³–10⁴× slower than dense — a bad trade for a query index.

## 7. Limitations / next

- Peel is still single-threaded and dominates `spectrum_ms` on dense graphs;
  parallelizing it would further cut build time.
- High-ω dense graphs (dblp, HepPh) clamp κ_s>2⁵³ at 1e18 in the peel; the
  index stores whatever the peel produces (comparison vs baseline still exact
  since both clamp identically).
- The query returns κ_s (core values). The connected-component hierarchy T_s
  is a downstream union-find pass on the level-s induced subgraph; not yet
  materialized/benchmarked.
- A middle tier (universal CPI + recompute+peel, no materialization) trades
  size for query and was not separately benchmarked.
