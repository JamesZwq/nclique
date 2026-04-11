# Region Tuple V2: Solving Open Problems

Reference: `docs/region_tuple_v1_complete.md` for V1 theory, code, and results.
Code: `src/NucleusDecomposition/NucleusCoreDecompositionRegion.cpp` (V1, frozen).

## Problems to Solve

### P1. Cascade for general s (s > r+1)

**Status**: V1 cascade only correct for s = r+1.

**Root cause**: V1 generates neighbor tuples by replacing 1 position.
For s = r+1: each (r+1)-clique has r+1 r-cliques, each differing by 1 vertex. ✓
For s > r+1: each s-clique has C(s,r) r-cliques, some differing by up to s-r vertices.

**Example**: r=3, s=5. A 5-clique {a,b,c,d,e} has C(5,3)=10 triangles.
Triangle {a,b,c} differs from {a,d,e} by 2 vertices (b,c → d,e).
V1 cascade would miss this connection.

**Solution: Bipartite Incidence between r-tuples and s-tuples**

Instead of computing cascade on-the-fly, enumerate BOTH r-tuples and s-tuples.
Both are compressed by the same region structure.

**s-tuple** σ = sorted multiset of s classes sharing a common clique.
Compression: C(|M|,s) s-cliques → C(ρ+s-1, s) s-tuples.

**Incidence**: For each s-tuple σ, list its r-sub-tuples (r-sub-multisets).
For each r-sub-tuple τ ⊂ σ: store ext(σ,τ) = Π_i C(|c_i| - j_i, m_i - j_i)
where j_i = multiplicity of c_i in τ, m_i = multiplicity in σ.

**Algorithm**:
```
1. Enumerate r-tuples (same as V1)
2. Enumerate s-tuples (same multiset enumeration, size s)
3. Build incidence: for each s-tuple σ, find its r-sub-tuples
4. support(τ) = Σ_{alive σ ∋ τ} ext(σ, τ)
5. Peeling:
   Pop min-support r-tuple τ
   For each s-tuple σ ∋ τ:
     If alive[σ]:
       alive[σ] = false
       For each τ' ≠ τ in σ:
         support[τ'] -= ext(σ, τ')
         Update bucket queue
```

**Why this is correct**: An s-tuple σ dies when ANY of its r-sub-tuples
is peeled. This exactly mirrors standard peeling where an s-clique becomes
incomplete when any r-clique is peeled. ext(σ,τ') correctly counts
s-clique instances.

**Complexity**: O(#s-tuples × avg-incident-r-tuples). Both compressed.
No individual r-clique or s-clique enumeration.

**Reduction formula**: ext(σ, τ') = Π_i C(|c_i| - j_i, m_i - j_i)
where j_i = class c_i's count in τ', m_i = count in σ.

### P2. Sparse graph performance

**Status**: V1 is slower than ST on com-dblp (0.2x) and web-Stanford (0.08x).

**Root cause**: Region size ≈ 1 → tuple count ≈ r-clique count → hash map
overhead dominates.

**Approach A: Auto-select**
Before tuple enumeration: compute total expected tuples = Σ C(ρ(M)+r-1, r).
If tuples > threshold (e.g., 2M): fall back to ST.
Simple, robust, no code change to tuple algorithm.

**Approach B: Hybrid**
For each maximal clique M: if ρ(M) ≈ |M| (no compression), use ST-style
per-r-clique peeling for M's r-cliques. Otherwise use tuples.
More complex but better worst-case.

**Approach C: Optimize data structures**
Replace hash map with sorted array or CSR for tuple storage.
Reduce per-tuple overhead. Helps but doesn't fix the fundamental
compression issue.

### P3. SDCT_MaxClique tagging failures — RESOLVED

**Investigation result**:
- email-Enron: **empty file** (0 bytes). Not an algorithm bug.
- email-Eu-core: tags 34K/42.7K maximal cliques (80%).
  BUT: Region Tuple gives **EXACT MATCH** with ST on 18 core levels.
  The incomplete tagging doesn't affect correctness.
- All other tested graphs: correct tagging and EXACT results.

**Conclusion**: P3 is NOT a blocking issue. No fix needed.

### P4. Skip unnecessary index builds — DONE

cliqueIndex and treeGraphV skipped when `regionOnly` is true.
Saves 5.3s on web-it-2004. Correctness verified.

### P5. SDCT build time dominance — DEFERRED

SDCT_MaxClique: 82s on web-it-2004 (96% of total time).
Most time is BK recursion, not tree storage.
Skipping sub-clique storage wouldn't help much.

**Real fix**: replace SDCT with a dedicated maximal clique enumerator
(e.g., quick-cliques in the repo). Region V2 only needs maximal clique
vertex sets, not CPI hold/pivot info. This is a bigger engineering task.

## Implementation Plan

New file: `src/NucleusDecomposition/NucleusCoreDecompositionRegionV2.cpp`
- Copy from V1 as starting point
- Address P1 (general s cascade) first — most important for paper
- Then P2 (auto-select) — practical impact
- P3-P5 are optimizations, lower priority

## Theoretical Work Needed

For P1: prove that the (s-r)-extension enumeration at the region level
gives correct cascade. Need to show:
1. Each s-clique is counted once
2. The cascade to each affected tuple is correct
3. The total cascade matches per-r-clique cascade

For the paper: the V1 theory (Theorems 1-4) is solid for s=r+1.
Extension to general s needs the same rigor.
