# r1Hier — Nucleus Spectrum Index: Handoff

**Last updated:** 2026-06-08
**Scope:** r=1 vertex s-clique-core *spectrum* index — "given any s in [2,ω],
return the s-clique core values κ_s(v) (and, eventually, the hierarchy T_s)".
Built on top of the SPIN★ (ST_V3) machinery in this repo.

> Data-safety note: `docs/` and `bench_plots/` are **gitignored**. The canonical
> copies of all findings docs + CSVs live here in `r1Hier/artifacts/` (git-tracked).
> All source code IS committed on `main` (see Commits below).

---

## 0. TL;DR status

| Piece | State |
|---|---|
| Math foundation (KK law, nesting, saturation) | verified, brute-force + on 9 real graphs |
| Saturation / compression measurement (`bench_spectrum_sat`) | done, 9 graphs |
| Build-once **universal CPI** (`SDCT_Augmented_NoTree_Universal`) | done, correct (0 mismatches incl. all-s) |
| Materialized index + query (`bench_spectrum_index`) | done: build/size/query/memory all measured |
| Hierarchy T_s materialization (components per level) | **NOT done** (next) |
| Peel parallelization | **NOT done** (single-thread; dominates build on dense graphs) |

**Headline results (single-thread, Darwin arm64):**
- Build-once universal CPI vs ω-rebuild: **1.6×–233×** faster (grows with ω).
- Dense index query: **µs latency**, **113×–88000×** faster than the *fair*
  build-once-tree baseline (recompute+peel), 869×–7.4M× vs full rebuild.
- Memory: dense index is **0.4–0.6× (down to 0.07×) the CPI tree** it was built
  from — materialize and drop the tree.
- Anchor+delta compression: 1.3–6× smaller than dense but query 10³–10⁴× slower
  → a bad trade for a query index.

---

## 1. The problem & the index

For each vertex v, the spectrum is its trajectory κ_2(v), κ_3(v), …, κ_{s*}(v),
where κ_s(v) = max k s.t. v is in a subgraph where every vertex has ≥ k
s-cliques, and s*(v) = max{s : κ_s(v) > 0}.

Two materialized index layouts:
- **dense-by-s**: per-level CSR of active `(vertex, κ)`. `query(s)` = read level-s
  slice. O(n_active_s), optimal. ~12 B/entry × Σ_s n_active_s.
- **anchor+delta**: per vertex `(s*, κ_{s*})` + sparse Kruskal-Katona deltas;
  `query(s)` reconstructs downward. Smaller, slow query.

---

## 2. Verified math foundation (do NOT re-derive wrong)

- **Shifted Containment Law (HOLDS, tight):** κ_s(v) ≥ ∂_s(κ_{s+1}(v)), where
  ∂_s = Kruskal-Katona **lower shadow** (s-cascade: write k = C(a_s,s)+
  C(a_{s-1},s-1)+…, shadow = C(a_s,s-1)+C(a_{s-1},s-2)+…). Verified law_viol=0
  on all 9 real graphs (exact regime).
  - WRONG earlier form g=C(d*,s-1) overclaims; only correct when k is exactly C(d,s).
- **Support nesting (HOLDS):** {v: κ_{s+1}>0} ⊆ {v: κ_s>0}; so active s-range is
  contiguous [2, s*]. nest_viol=0 on all graphs.
- **Non-monotone in s:** κ_s NOT ≥ κ_{s+1}. K_m vertex: κ_s = C(m-1,s-1), e.g.
  K_6 → [5,10,10,5,1], K_5 → [4,6,4,1] (these are fully KK-saturated).
- **δ_s = κ_s − ∂_s(κ_{s+1}) ≥ 0**; saturated vertex = all δ=0 (trajectory =
  anchor alone). `kk_shadow` self-test: `bench_spectrum_index --selftest`.
- **3 RETRACTED overclaims** (don't state as theorems): (1) Ω(nω) space lower
  bound is not info-theoretic; (2) "must do ω independent peels" unproven;
  (3) cross-s pair query is a TIE not a win. See memory file.

---

## 3. Critical engineering findings

1. **Per-s SDCT is NOT reusable across s.** `SDCT_Augmented_NoTree(g, max_k=s)`
   emits a leaf only when maximal-clique cSize ≥ s (`bkRecurse_NoTree`, line 51),
   with per-s (keepV,dropV). So "build once at ω, recompute all s" is INVALID
   with that function.
2. **Fix = universal SDCT.** `SDCT_Augmented_NoTree_Universal` (added) emits
   EVERY maximal-clique leaf once (no max_k prune, no cSize≥max_k filter). One
   traversal serves all s via the Pivoter invariant: leaf (H,P) contributes
   C(|P|, s−|H|) s-cliques; per-vertex support at level s is C(|P|,s−|H|) for
   v∈H and C(|P|−1,s−|H|−1) for v∈P. Then reuse the validated `ST_V3_Peel`.
3. **Numerical clamp.** On high-ω dense graphs (dblp ω=114, HepPh ω=239) κ_s
   exceeds 2⁵³ and the peel clamps support at 1e18 (`kBucketKeyClamp`). The
   high-s spectrum is not double-exact; flagged `fuzzy`, excluded from trustworthy
   stats. NOT counted as law violations. (On arm64 `long double == double`.)
4. **Fair query baselines (3 tiers):** full per-s rebuild (R builds) < build
   universal tree once + recompute+peel per query (shared) < dense materialized.
   Report dense vs *shared* as the honest headline.

---

## 4. Results (full data in `r1Hier/artifacts/spectrum_index.csv`, `spectrum_sat.csv`)

### Build: universal build-once vs ω-rebuild
| graph | ω | index build ms | ω-rebuild ms | speedup |
|---|---|---|---|---|
| amazon0302 | 7 | 452 | 765 | 1.7× |
| ca-CondMat | 26 | 59 | 286 | 4.9× |
| ca-AstroPh | 57 | 360 | 4742 | 13.2× |
| com-dblp | 114 | 931 | 29863 | 32.1× |
| ca-HepPh | 239 | 1290 | 211379 | 163.8× |

### Query latency (s=2, the least favorable; high-s reaches 10⁴–10⁶×)
| graph | dense µs | shared ms | rebuild ms | dense/shared | dense/rebuild |
|---|---|---|---|---|---|
| ca-CondMat | 21.3 | 4.48 | 24.7 | 211× | 1160× |
| ca-AstroPh | 20.7 | 10.55 | 160.6 | 510× | 7758× |
| com-dblp | 296.5 | 33.6 | 379.5 | 113× | 1280× |
| ca-HepPh | 13.5 | 4.42 | 574.7 | 328× | 42621× |

### Memory: index vs CPI "tree" substrate
| graph | CPI tree MB | dense MB | a+d MB | dense/tree | a+d/tree |
|---|---|---|---|---|---|
| ca-CondMat | 2.59 | 1.57 | 0.49 | 0.61× | 0.19× |
| ca-AstroPh | 8.14 | 3.15 | 0.52 | 0.39× | 0.06× |
| cit-HepPh | 46.24 | 3.15 | 1.82 | 0.07× | 0.04× |
| com-dblp | 26.36 | 25.17 | 6.50 | 0.95× | 0.25× |

### Saturation / compression (anchor+delta, from `bench_spectrum_sat`)
- Saturation fraction 13%–89% (collaboration nets high, citation low).
- Compression comp = avgTraj/(2 + avg δ≠0), capped at avgTraj/2; range
  **0.90× (cit-HepPh, EXPANDS) to 3.67× (ca-AstroPh)**. Modest.

---

## 5. Code & artifacts

**Source (committed on `main`):**
- `src/SDCT_Augmented.h` / `.inl` — `SDCT_Augmented_NoTree_Universal` (build-once).
- `bench_spectrum_sat.cpp` — ω + saturation/compression + KK-law validation.
- `bench_spectrum_index.cpp` — materialized index: build/size/query/memory.
  Backends: dense-by-s, anchor+delta. `IDX_BUILD/IDX_SIZE/IDX_MEM/IDX_VERIFY/IDX_QUERY` lines.
- `run_spectrum_sat.sh` — batch driver → CSV.

**Reused (do not re-implement):** `NCliqueVertexCoreDecomposition_ST_V3_Build/_Peel`
(`src/NucleusDecomposition/NCliqueVertexCoreDecompositionST_V3.cpp`), the trusted
SPIN★ path (bit-identical to `degeneracy_cliques` under `PIVOTER_RUN_ST_V3`).

**Artifacts (git-tracked copies here; originals in gitignored docs/bench_plots):**
- `r1Hier/artifacts/spectrum_index_findings.md`, `spectrum_saturation_findings.md`
- `r1Hier/artifacts/spectrum_index.csv`, `spectrum_sat.csv`
- `r1Hier/artifacts/nucleus_spectrum_animation.html`

**Memory file:** `~/.claude/projects/.../memory/project_spectrum_index_core.md`
(verified math, results, retracted overclaims).

---

## 6. Build & reproduce

```bash
cd /Users/zhangwenqian/UNSW/pivoter
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 12 --target bench_spectrum_index bench_spectrum_sat

# self-test the KK shadow
./build/bin/bench_spectrum_index --selftest

# full index experiment on one graph (build/size/mem/verify/query)
./build/bin/bench_spectrum_index graphs/ca-CondMat.edges
# exhaustive all-s correctness (no perf):
./build/bin/bench_spectrum_index graphs/ca-GrQc.edges \
    --queries $(python3 -c "print(','.join(map(str,range(2,45))))") --no-baseline

# saturation/compression sweep:
./build/bin/bench_spectrum_sat graphs/ca-AstroPh.edges
```
- Build cap: `-j 12` (hard repo rule).
- Test graphs in `graphs/*.edges`. Avoid `facebook_combined.edges` (too dense).
- Correctness gate ALWAYS before perf claims: `IDX_VERIFY ... mismatches=0`.

---

## 7. Caveats / honest limits

- The dense index answers **only κ_s**; the CPI tree can also recompute any s and
  support other queries (clique counts, etc.). "Drop the tree" assumes κ_s-only.
- Hierarchy **T_s (connected components per level) is NOT materialized** yet —
  query currently returns core values, not the nested component tree.
- Peel is single-thread and dominates `spectrum_ms` on dense graphs.
- Numerical clamp on κ_s > 2⁵³ (high-ω dense graphs); high-s spectrum approximate
  there (both index and baseline clamp identically, so comparisons stay exact).
- Compression (anchor+delta) is modest (≤~3.7×, can expand) — weak standalone story.

See `r1Hier/TODO.md` for prioritized next steps.
