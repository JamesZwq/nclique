# Region Tuple (r,s)-Nucleus Decomposition — V1 Complete Record

## 1. Theory

### 1.1 Definitions

- Graph G=(V,E), maximal clique set M = {M_1, ..., M_m}
- **Profile**: prof(v) = {M ∈ M : v ∈ M}
- **Region (Overlap Class)**: R = maximal set of vertices with identical profile
  - p ≤ n regions, disjoint partition of V
- **ρ(M)** = number of regions in clique M, **α = |M|/ρ** = avg region size
- **Region r-tuple** τ = sorted multiset of r regions sharing a common clique
- **Multiplicity**: mult(τ) = Π_R C(|R|, τ_R)

### 1.2 Proven Theorems

**Lemma 1 (Completeness)**: Every r-clique is represented by exactly one tuple.

**Lemma 2 (Containment)**: For T ⊆ R (same-region r-clique):
{M : T ⊆ M} = prof_R. Proof: 3 lines via profile identity.

**Theorem 1 (Support Equality)**: All r-cliques in the same tuple have
identical support for any s ≥ r+1. Proof: containing cliques = ∩ prof_{R_i},
depends only on tuple, not specific vertices.

**Theorem 2 (Compression)**: compression ≈ α^r per clique.
Grows exponentially with r.

**Theorem 3 (Tuple count)**: C(ρ+r-1, r) per clique = O(r^{ρ-1}).
Polynomial in r, vs factorial for r-cliques.

**Theorem 4 (Equivalence)**: Tuple peeling gives EXACT same core numbers
as per-r-clique peeling (for s=r+1, verified experimentally).

### 1.3 Support Formula

For s = r+1:  support(τ) = |∪_{M ∈ C_τ} alive_M| - r

For general s:  support(τ) = Σ_{∅≠S⊆C_τ} (-1)^{|S|+1} C(|∩_{M∈S} alive_M| - r, s-r)

Where C_τ = ∩_j prof_{R_{i_j}} = common cliques of the tuple.

### 1.4 Cascade (s = r+1)

When tuple τ is peeled, for neighbor class D and position k:
- Neighbor tuple τ' = sort(τ with position k replaced by D)
- Reduction = |τ[k]| - |{j≠k : τ[j] = τ[k]}|
- Check: other r-1 neighbor tuples (position j≠k replaced by D) must be alive
- Core assignment: coreLevel = max(coreLevel, currentBucketLevel)  ← NON-DECREASING

## 2. Implementation

### 2.1 Files

- `src/NucleusDecomposition/NucleusCoreDecompositionRegion.cpp` — main implementation
- `src/NucleusDecomposition/NucleusCoreDecompositionRegion_backup.cpp` — pre-bugfix backup
- `src/SDCT_MaxClique.hpp` — SDCT with X tracking for maximal clique tagging
- `src/degeneracy_cliques.cpp` — dispatch via PIVOTER_RUN_REGION env var (r≥3)

### 2.2 Algorithm Steps

1. **Build maximal-clique regions**: SDCT_MaxClique tags maximal vs sub-clique paths.
   Merge sub-cliques into parent maximal cliques. O(L).

2. **Build overlap classes**: Group vertices by maximal-clique membership profile.
   classOf[v] = class ID. O(n).

3. **Enumerate region r-tuples**: For each maximal clique M, recursive multiset
   enumeration of r classes from regions(M). O(Σ C(ρ(M)+r-1, r)).

4. **Compute support**: For s=r+1: union of common clique vertices - r.
   For general s: inclusion-exclusion (2^|C_τ| terms per tuple).

5. **Cascade peeling**: Bucket queue, peel minimum support tuple.
   For each neighbor class D × position k: compute neighbor tuple, check completeness,
   apply reduction. coreLevel only increases.

### 2.3 Bug Found and Fixed

**Bug**: `currentLevel` could decrease during cascade → core values assigned
at wrong (lower) levels.

**Fix**: Added `coreLevel = max(coreLevel, currentLevel)`. Core values are
now non-decreasing, matching the standard peeling invariant.

**Verification**: Exact match on all tested graphs after fix.

## 3. Experimental Results

### 3.1 Correctness (s = r+1): ALL EXACT

| Graph | Vertices | r | s | Core Levels | ST vs Tuple |
|-------|----------|---|---|-------------|-------------|
| debug_small | 8 | 3 | 4 | 2 | EXACT |
| dblp-core30 | 1,206 | 3 | 4 | 29 | EXACT |
| dblp-core30 | 1,206 | 4 | 5 | 28 | EXACT |
| com-dblp | 317,080 | 3 | 4 | 45 | EXACT |
| web-Stanford | 281,903 | 3 | 4 | 58 | EXACT |

Verification method: line-by-line diff of all core levels and counts.

### 3.2 Correctness (s > r+1): MISMATCH

| Graph | r | s | Issue |
|-------|---|---|-------|
| dblp-core30 | 3 | 5 | Diff at core=171 (cascade incomplete) |

Root cause: cascade only generates 1-position replacements (correct for s=r+1).
For s>r+1: s-cliques have C(s,r) r-cliques, some differing by >1 position.

### 3.3 Performance

| Graph | r | s | ST time | Tuple time | Tuples | Speedup |
|-------|---|---|---------|-----------|--------|---------|
| dblp-core30 | 3 | 4 | 380ms | 4ms | 1,706 | **95x** |
| dblp-core30 | 4 | 5 | 13.3s | 17ms | 4,434 | **784x** |
| dblp-core30 | 5 | 6 | 566s | 33ms | 7,211 | **17,000x** |
| com-dblp | 3 | 4 | 809ms | 3,613ms | 934K | **0.2x** (slower) |
| web-Stanford | 3 | 4 | 17.7s | 211s | 6.3M | **0.08x** (slower) |
| web-it-2004 | 3 | 4 | >1h TO | 3.2s | 788K | **>1000x** |

### 3.4 When Region Tuple Wins vs Loses

Wins (α >> 1): dblp-core30 (α≈12), web-it-2004 (α≈9)
Loses (α ≈ 1): com-dblp, web-Stanford (most regions size 1, tuple ≈ r-clique count)

## 4. Open Problems for V2

### P1. Cascade for general s (s > r+1)

Current cascade replaces 1 position → generates r neighbor tuples.
For s>r+1: need to replace up to s-r positions simultaneously.
Each s-clique has C(s,r) r-cliques, some differing by >1 vertex.

### P2. Sparse graph performance

On sparse graphs: region size ≈ 1 → no compression → slower than ST
(hash map overhead). Options:
- Auto-select: if avg region size < threshold → fall back to ST
- Optimize hash map → flat array when tuple count is small
- Hybrid: use tuples for large cliques, direct r-cliques for small cliques

### P3. SDCT_MaxClique tagging failures

email-Enron: tags 0 maximal cliques. email-Eu-core: tags 34K vs 42.7K (NetworkX).
The X propagation in BK has bugs for some graph structures.

### P4. Clique Index build overhead

The cliqueIndex mapList build takes significant time (5.3s on web-it-2004).
Region Tuple doesn't need this — could skip it.

### P5. SDCT build time dominance

SDCT_MaxClique: 87s on web-it-2004 (vs 3.2s for tuple decomposition).
For large r: SDCT builds for s=r+1 become very expensive (memory + time).
Could extract maximal clique info without storing full tree.
