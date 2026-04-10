# CPI-Path-Based (r,s)-Nucleus Decomposition: Theory v2

## 0. Key Insight

CPI paths (SDCT leaves) partition ALL s-cliques. Number of paths is a
GRAPH PROPERTY, independent of r. If we can do nucleus decomposition
at the PATH level: complexity is independent of r.

## 1. CPI Path Structure (Review)

Each CPI path P = (V_h(P), V_p(P)):
- V_h = hold vertices (|V_h| = h)
- V_p = pivot vertices (|V_p| = p)
- Path size: |P| = h + p
- P is a CLIQUE (all vertices mutually adjacent)
- P generates exactly C(p, s-h) distinct s-cliques
- Each s-clique uses ALL h hold vertices + (s-h) of the p pivots
- **Key**: s-cliques are PARTITIONED across paths (each s-clique in exactly one path)

## 2. Region := CPI Path

**Definition**: Region R_i = CPI path P_i.

Properties:
- Regions OVERLAP (a vertex can be in multiple regions)
- Every region is a clique
- Number of regions = number of CPI paths (independent of r)
- web-it-2004: 424K regions for ANY r

## 3. R-Clique Support via Paths

For r-clique T: support(T) = number of s-cliques containing T.

Since s-cliques are partitioned by paths:

  support(T) = Σ_{P : T ⊆ P} C(p_P - b_P(T), s - h_P - b_P(T))

where b_P(T) = |T ∩ V_p(P)| = number of T's vertices that are pivots in P.

This is the standard Lemma 4 formula.

## 4. Path-Level Decomposition: Ideas

### 4.1 Approach: Per-Path Contribution Tracking

For each path P, maintain counters:
- alive_h(P) = number of alive hold vertices
- alive_p(P) = number of alive pivot vertices
- s_cliques(P) = C(alive_p, s - alive_h) = number of alive s-cliques from P

When a vertex v is removed from the graph:
- For each path P containing v:
  - If v is hold in P: alive_h(P) -= 1
  - If v is pivot in P: alive_p(P) -= 1
  - Recompute s_cliques(P)

This gives us the total number of alive s-cliques at any point.

### 4.2 But What Gets Peeled?

Standard (r,s)-nucleus: peel r-cliques by support.
Path-level: what do we peel?

Options:
A) Peel VERTICES. When vertex v is peeled: remove all r-cliques involving v.
   - Pro: vertex count is independent of r
   - Con: not the standard (r,s)-nucleus (which peels r-cliques, not vertices)

B) Peel PATHS. When path P is peeled: remove all its s-cliques.
   - Pro: path count is independent of r
   - Con: different from standard nucleus

C) Peel R-CLIQUE GROUPS. Group r-cliques by (containing paths, subNumPivot vector).
   - Pro: exact standard nucleus
   - Con: number of groups might depend on r

### 4.3 Vertex-Level Peeling with Path Counters

This is analogous to the R=1 ST approach, generalized to r ≥ 2:

For each vertex v: 
  vertex_support(v) = Σ_{r-cliques T involving v} support(T)

This is a SUM over all r-cliques involving v. Can we compute it from path counters?

vertex_support(v) = Σ_{P containing v} Σ_{T ⊆ P, v ∈ T} C(p_P - b_P(T), s - h_P - b_P(T))

The inner sum = sum over r-cliques in P that include v, weighted by their per-P contribution.

For vertex v in path P (v is hold with index in H):
  Σ_{T ∋ v, T ⊆ P} C(p-b, s-h-b) = Σ_b C(h-1, r-1-b) · C(p, b) · C(p-b, s-h-b)
                                     (choosing r-1-b more holds and b pivots for T)

For vertex v as pivot:
  Σ_{T ∋ v, T ⊆ P} C(p-b, s-h-b) = Σ_b C(h, r-b) · C(p-1, b-1) · C(p-b, s-h-b)
                                     (choosing r-b holds and b-1 more pivots)

These are computable from (h, p, r, s) in O(r) time per path per vertex.

But: vertex_support ≠ the support of any single r-clique. It's a SUM.
Peeling by vertex_support gives a DIFFERENT decomposition from r-clique peeling.

### 4.4 TODO: What Exactly Is the Right Path-Level Metric?

OPEN QUESTION: What path-level or vertex-level metric gives a decomposition
that is:
1. Equivalent to standard (r,s)-nucleus, OR
2. A provably related structural decomposition, AND
3. Computable in time independent of r?

## 5. What We Know For Sure

### 5.1 CPI Enables r-Independent Computation of TOTALS

- Total s-cliques = Σ_P C(p_P, s - h_P). O(#paths).
- Total support sum = C(s,r) · total_s_cliques. O(#paths).
- Per-vertex total support = sum of per-path formulas. O(#paths × r).

All independent of the number of r-cliques.

### 5.2 CPI Does NOT Directly Enable r-Independent PEELING

Standard peeling requires per-r-clique supports. An r-clique's support
depends on its specific path profile. Different r-cliques from the same
path can have different total supports.

### 5.3 Region Tuple Gives Exact Decomposition

Region tuples (Section 4 of theory v1) give EXACT (r,s)-nucleus.
Number of tuples << number of r-cliques on dense graphs.
But tuple count depends on r: O(Σ_M C(ρ(M), r)).

## 6. Possible Resolution: R=1-Style Counter Approach

For R=1: each vertex has support = Σ_P w(v, P). When a vertex is peeled:
per-path counters update, and support changes by a delta.

For R ≥ 2: can we define a "quotient class" that avoids r-clique enumeration?

Within each path P: r-cliques with the same subNumPivot b form a quotient class.
Number of classes per path: O(r+1).
Total classes: O(#paths × r).

An r-clique's total support = Σ_P contribution(P, b_P).
Two r-cliques in the SAME quotient class in EVERY path have the same support.

The number of distinct "global quotient classes" (same class in every path)
is bounded by O(product of (r+1) over all containing paths) per r-clique.

For most r-cliques in 1-2 paths: O(r+1) classes. Manageable.
For dense-core r-cliques in many paths: could be large.

## 7. Summary

| Approach | Exact? | r-independent? | Entity count |
|----------|--------|----------------|-------------|
| Standard (r-clique) | Yes | No | Σ C(|M|, r) |
| Overlap class (region-only) | No (overestimates) | Yes | ≤ n |
| Region tuple (v1) | Yes | No | Σ C(ρ(M), r) |
| CPI path (v2) | ? | ? | #paths (fixed) |

The CPI path approach is promising but needs a precise formulation of
what metric to use for peeling.
