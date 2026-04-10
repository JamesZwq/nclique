# Region-Based Nucleus Decomposition

## Core Idea

Instead of tracking individual r-cliques (338M on web-it-2004),
decompose the vertex set into **overlap regions** defined by maximal clique membership.

All r-cliques within the same region have the same support and core value.
The complexity depends on the number of regions, not the number of r-cliques.

## Definitions

### Vertex Profile

For each vertex v, define:

```
profile(v) = { (pathId, isPivot) : v is in path pathId }
```

Two vertices with the same profile appear in exactly the same set of paths
with the same hold/pivot status.

### Region

A **region** is a maximal set of vertices sharing the same profile.

Properties:
- Every r-clique with all vertices from the same region has the same support
- Region size = k vertices → C(k, r) r-cliques, all with the same core value
- Number of regions ≤ number of vertices (typically much smaller)
- **Does not depend on r**

### Region Support

For a region R with profile P = {(path_1, isPivot_1), ..., (path_m, isPivot_m)}:

An r-clique from R with subNumPivot = q in path path_j has
local support = C(pivotC_j - q, s - holdC_j - q) (Lemma 4 from paper).

Total support = Σ over all paths in the profile of the local support.

Since all vertices in R have the same hold/pivot status in each path,
q is determined by the region and the number of pivot vertices chosen.

## Why This Works on web-it-2004

- 509K vertices, 424K paths
- Dense 430-core: ~430 vertices in ~400+ paths each
- Most core vertices share the same profile → few large regions
- Peripheral vertices have unique profiles → many tiny regions (size 1-3)
- Total regions: O(thousands), not O(338M r-cliques)

## Complexity

- Region computation: O(vertices × avg_paths_per_vertex) — one-time
- Support computation: O(regions × paths_per_region) — no r-clique enumeration
- Peeling: O(regions) — each region is peeled as a unit
- Triangle count per region: C(region_size, r) — just arithmetic, no enumeration

## Step 1 Result: Pivot-Profile Regions

| Graph | Vertices | Regions | Max region | Time |
|---|---|---|---|---|
| email-Eu-core | 806 | 532 | 274 | 1 ms |
| dblp-core30 | 1,206 | 1,142 | 32 | 0 ms |
| web-it-2004 | 437,936 | 360,807 | 31,033 | 125 ms |

## Step 2 Result: Same-Region + Cross-Region Triangle Count

| Graph | Same-region | Cross-region | Total | V20 total | Coverage |
|---|---|---|---|---|---|
| dblp-core30 | 1,747 | 673,957 | 675,704 | 686,193 | 98.5% |
| web-it-2004 | 997,970 | 336,065,299 | 337,063,269 | 338,476,979 | 99.6% |

Missing ~0.4-1.5%: triangles where one vertex is a HOLD in a path but NOT
the earliest vertex. These need separate handling in Step 3.

## Steps 4+5 Result: Region-Triples + Static Peeling

| Graph | Region-triples | Triangles | Max core | Time |
|---|---|---|---|---|
| dblp-core30 | 645K | 675K | 159 | 59 ms |
| web-it-2004 | **327M** | 337M | 430 | **92.6s** |

### Critical finding: region-triples ≈ individual triangles

On web-it-2004: 327M region-triples vs 338M triangles. Almost no compression!

Root cause: 93% of regions have size=1 (only one vertex). Cross-region triples
between two size=1 regions are just individual triangles. The strict pivot-profile
requirement makes most regions trivially small.

### Why pivot-profile regions are too strict

Two vertices with pivotProfiles differing by even ONE path are in different regions.
On web-it-2004 with 424K paths, most vertices differ by at least one path → size=1.

The 31K-vertex region (pivotPaths=0) is the exception: vertices that don't appear
as pivot in ANY path. These are "leaf-only" vertices.

### Next direction: coarser grouping

Need a grouping that creates LARGER regions. Options:

1. **Group by pivotC only**: vertices with same pivotC value → same region.
   Simple but loses too much information (many different vertices have same pivotC).

2. **Group by "core membership"**: vertices that will be peeled at the same core level
   form a group. But this is circular (core level is what we're computing).

3. **Group by hold-path pivotC + pivot-path count**: coarser than exact pivot profile.

4. **Hierarchical regions**: start with strict profiles, merge regions that share
   most paths. Merge until region sizes are large enough.

5. **CPI subtree regions**: group by which subtree of the degeneracy tree the
   vertex belongs to. Vertices in the same subtree share most paths.

## Implementation Progress

1. ~~Build pivot-profile regions~~ ✓ (but too strict)
2. ~~Compute same-region + cross-region support~~ ✓ 
3. ~~Handle non-keepC=1 paths~~ ✓ (identified 1-2% missing)
4. ~~Build region-triples + static peeling~~ ✓ (but 327M triples = no compression)
5. Cascade peeling — blocked by region size issue
6. **Need: coarser region definition for real compression**

## Maximal Clique Analysis (NetworkX ground truth)

| Graph | Maximal cliques | CPI paths | Sub-clique % | Max clique size |
|---|---|---|---|---|
| dblp-core30 | **87** | 1,230 | **93%** | 114 |
| email-Eu-core | 42,709 | 43,220 | 1.2% | 18 |
| web-it-2004 | 381,110 | 423,886 | 10% | 432 |

### Key finding

On dblp-core30: 87 maximal cliques vs 1,230 CPI paths = **14x compression**.
Largest maximal clique has 114 vertices → C(114, 3) = 242,556 r-cliques in ONE region.

CPI produces sub-cliques because it removes BK's X (exclusion set) for counting.
Restoring X tracking (without changing path construction) tags each path as maximal or sub.

### SDCT_MaxClique implementation

Added X tracking to BK recursion in `src/SDCT_MaxClique.hpp`:
- X initialized as earlier neighbors at top level
- X ∩ N(v) propagated through recursion
- Path tagged maximal iff P=∅ AND X=∅

Current bug: tags 68 maximal on dblp-core30 (should be 87). X propagation
has a false-positive issue (19 true maximal cliques incorrectly tagged as sub).

## CRITICAL: Support formula verification FAILED

Tested on email-Eu-core r=3 s=4 (104,730 triangles):

- **99.4% of triangles have NO keepC=1 "earliest" hold-path** containing the other two
- Of the 0.6% found: **33% have WRONG support** (region formula >> ground truth)
- Example: gt=44, region=360 (8x overestimate)

### Root cause

The assumption "every triangle {a,b,c} with a earliest has a's hold-path containing
b,c as pivots" is **WRONG**. CPI paths are NOT maximal cliques — they're BK search
tree paths. A CPI path can have holdC > 1, and many triangles only appear in
holdC > 1 paths.

The formula `support = contribP2(a) + n_allPivotPaths` is also wrong because:
1. It assumes the triangle has exactly ONE "hold-path" contribution — but it can
   have contributions from MULTIPLE paths with different hold/pivot structures
2. The "contribP2" formula assumes keepC=1 — but paths can have keepC > 1

### Root cause (from re-reading the paper)

The Build-CPI algorithm (Algorithm 3) constructs paths via BK search.
At each BK branch: the BK pivot → CPI pivot, non-pivot branch → CPI hold.
So paths can have holdC = 1, 2, 3, ... (one hold per BK branch level).

On email-Eu-core: most paths have holdC > 1. 99.4% of triangles are ONLY
in holdC > 1 paths. The "keepC=1 earliest vertex" model covers < 1%.

### Fundamental limitation of vertex-profile regions

Support depends on PATH-LEVEL hold/pivot assignment, not vertex-level properties.
The same triangle in different paths has different hold/pivot splits → different
per-path contributions. No vertex-level grouping can capture this.

### What remains correct

- **Per-path quotient classes** (Lemma 4, per subNumPivot): CORRECT within each path
- **Per-leaf accumulation** (Σ over paths of C(p-b, s-h-b)): CORRECT for total support
- **DPCM locality** (only affected paths change): CORRECT for cascade

### What's wrong

- "Vertex pivot profile determines support": WRONG (hold/pivot varies per path)
- "keepC=1 is the common case": WRONG (most paths have keepC > 1 on email-Eu-core)
- "contribP2(a) + n_allPivotPaths" as support formula: WRONG (misses holdC > 1 contributions)

### Next direction

The region approach cannot work at the vertex level. It must work at the
PATH-CLIQUE level: group (path, r-clique-class) pairs, not vertices.
This is exactly what the original quotient approach does.

The fundamental challenge remains: how to compute TOTAL support (sum across paths)
without enumerating individual r-cliques. The per-path quotient gives excellent
within-path compression (34,055x on web-it-2004), but cross-path aggregation
still requires per-r-clique tracking (or equivalent).
