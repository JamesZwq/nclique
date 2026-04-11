# Region Tuple V2: Complete Record

## Algorithm Overview

Region Tuple decomposes (r,s)-nucleus without enumerating individual r-cliques.
Group structurally equivalent r-cliques into "region tuples" and peel tuples.

### Step 1: Enumerate Maximal Cliques
MaxCliqEnum: degeneracy-ordered BK with P∪X pivoting.
Output: list of maximal cliques (vertex sets).

### Step 2: Build Overlap Classes (Regions)
prof(v) = {maximal cliques containing v}.
Region = vertices with identical prof. Region count p ≤ n.

### Step 3: Enumerate Region r-Tuples and s-Tuples
r-tuple = multiset of r regions sharing a common clique.
s-tuple = multiset of s regions sharing a common clique.
Build bipartite incidence: each s-tuple lists its r-sub-tuples.

### Step 4: Compute Initial Support
support(τ) = Σ_{alive σ ∋ τ} ext(σ, τ)
ext(σ, τ) = Π_i C(|c_i| - j_i, m_i - j_i)

### Step 5: Cascade Peeling
Bucket queue. Pop min-support r-tuple. Mark incident s-tuples dead.
Update neighbor r-tuples' support. coreLevel only increases.

## Theorems

### Theorem 1 (Containment Lemma)
For T ⊆ R: {M : T ⊆ M} = prof_R.
Proof: 3 lines. Any M containing T must be in prof of T's vertices = prof_R.

### Theorem 2 (Support Equality)
All r-cliques in the same tuple have identical support for any s.
Proof: containing cliques = ∩ prof_{R_i}, depends only on tuple.

### Theorem 3 (Completeness)
Every r-clique belongs to exactly one tuple.

### Theorem 4 (Compression)
Per clique M: compression ≈ α^r where α = |M|/ρ(M).

### Theorem 5 (Equivalence)
Tuple peeling gives EXACT same core numbers as per-r-clique peeling.

## Correctness Verification

11 test cases, ALL EXACT (line-by-line diff of every core level and count):

| Graph | r | s | Core Levels | Match |
|-------|---|---|-------------|-------|
| dblp-core30 | 3 | 4 | 29 | EXACT |
| dblp-core30 | 3 | 5 | 28 | EXACT |
| dblp-core30 | 3 | 6 | 27 | EXACT |
| dblp-core30 | 4 | 5 | 28 | EXACT |
| dblp-core30 | 4 | 6 | 27 | EXACT |
| dblp-core30 | 5 | 6 | 27 | EXACT |
| com-dblp | 3 | 4 | 45 | EXACT |
| com-dblp | 3 | 5 | 47 | EXACT |
| email-Eu-core | 3 | 4 | 18 | EXACT |
| email-Eu-core | 3 | 5 | 115 | EXACT |
| email-Eu-core | 4 | 5 | 16 | EXACT |

web-it-2004: ST/REF timeout, verified support sum + max core only.

## Performance

| Graph | r | s | REF | ST | V2 | V2 vs REF |
|-------|---|---|-----|-----|-----|-----------|
| dblp-core30 | 3 | 4 | 516ms | 411ms | 11ms | **47x** |
| dblp-core30 | 4 | 5 | 12.1s | 14.1s | 18ms | **672x** |
| dblp-core30 | 5 | 6 | >10min | >10min | 30ms | **>20,000x** |
| com-dblp | 3 | 4 | 1.27s | 870ms | 4.55s | 0.3x (slower) |
| web-it-2004 | 3 | 4 | >1h | >1h | 3.8s | **>1000x** |

MaxCliqEnum vs SDCT_MaxClique (web-it-2004): 8.3s vs 82s = 10x faster.
MaxCliqEnum finds 381K maximal cliques (complete), SDCT found 85K (22%).

## Implementation Files

| File | Description |
|------|-------------|
| `src/NucleusDecomposition/NucleusCoreDecompositionRegionV2.cpp` | V2: general (r,s), bipartite incidence |
| `src/NucleusDecomposition/NucleusCoreDecompositionRegion.cpp` | V1: s=r+1 only (frozen) |
| `src/MaxCliqEnum.hpp` | Lightweight maximal clique enumerator |
| `src/SDCT_MaxClique.hpp` | SDCT with X tracking (legacy, replaced by MaxCliqEnum) |
| `src/degeneracy_cliques.cpp` | Dispatch: PIVOTER_RUN_REGION_V2=1 |

## Usage

```bash
PIVOTER_RUN_REGION_V2=1 ./build/bin/degeneracy_cliques <graph.edges> <r> <s>
```

## Key Design Decisions

1. **Bipartite r-tuple ↔ s-tuple incidence** (not V1's per-position cascade):
   Handles any s correctly. s-tuple alive/dead naturally tracks s-clique completeness.

2. **MaxCliqEnum instead of SDCT**: P∪X pivoting + sub-clique pruning.
   10x faster, finds all maximal cliques (not just 22%).

3. **coreLevel = max(coreLevel, currentLevel)**: Standard peeling invariant.
   Bug found in initial implementation where currentLevel could decrease.

4. **Filter paths by size ≥ s** (not ≥ r): Matches ST baseline behavior.
   Small maximal cliques (size < s) contribute zero-support r-cliques.

## Limitations

1. Slower on sparse graphs (region size ≈ 1 → no compression).
2. s-tuple count can be large for big s values.
3. MaxCliqEnum BK implementation not as optimized as production enumerators.
