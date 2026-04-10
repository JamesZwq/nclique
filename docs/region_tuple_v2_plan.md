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

**Approach**: For general s, when peeling tuple τ:
- For each s-clique S containing an r-clique from τ:
  - S is determined by C_τ (common cliques) and choosing s-r more vertices
  - The C(s,r)-1 other r-cliques of S are in various tuples
  - Each gets -1 to support IF all OTHER r-cliques of S are alive

The challenge: enumerating the (s-r)-vertex extensions and mapping them back
to tuples, without enumerating individual r-cliques.

**Key insight**: Within a maximal clique M, an s-clique S uses vertices from
various regions. The r-cliques of S = all C(s,r) r-subsets. Each r-subset
maps to a tuple. The mapping depends on which regions the s-r extension
vertices come from.

For a tuple τ in clique M: the s-cliques containing τ's r-cliques use
s-r additional vertices from M. These vertices come from regions in
regions(M). For each multiset of (s-r) regions (the extension):
the s-clique's other r-cliques are determined.

This is a region-level enumeration of (s-r)-extensions, NOT per-r-clique.
Total: O(C(ρ(M)+s-r-1, s-r)) extensions per tuple, independent of |M|.

TODO: formalize and implement.

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

### P4. Skip unnecessary index builds

**Status**: cliqueIndex mapList build takes 5.3s on web-it-2004 but isn't
needed for Region Tuple.

**Approach**: When PIVOTER_RUN_REGION is set, skip the clique index build
entirely. Only build the SDCT tree and extract maximal clique tags.

### P5. SDCT build time dominance

**Status**: SDCT_MaxClique takes 87s on web-it-2004 (vs 3.2s for tuple decomp).

**Approach**: The Region Tuple algorithm only needs:
1. The set of maximal cliques (vertex sets)
2. Not the full SDCT tree with all sub-clique paths

Could build a "maximal-clique-only" SDCT that stops building sub-paths
once maximal status is determined. Or use an external maximal clique
enumerator (which is typically faster than full SDCT).

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
