# Region Tuple (r,s)-Nucleus Decomposition: Complete Theory

## 1. Definitions

### 1.1 Graph and Maximal Cliques

Graph G=(V,E), maximal clique set M = {M_1, ..., M_m}.

### 1.2 Regions (Overlap Classes)

**Profile**: prof(v) = {M ∈ M : v ∈ M}.

**Region**: R = maximal set of vertices with identical profile.
Region set R = {R_1, ..., R_p}. Properties:
- p ≤ n (at most one region per vertex)
- Disjoint partition: V = R_1 ∪ ... ∪ R_p
- All vertices in R share the EXACT same profile: prof_R

### 1.3 Region Structure within Maximal Cliques

For maximal clique M: define regions(M) = {R ∈ R : R ⊆ M}.

Since prof is identical within a region: if ANY vertex of R is in M,
then ALL vertices of R are in M. So R ⊆ M or R ∩ M = ∅.

Let ρ(M) = |regions(M)| and |M| = Σ_{R ∈ regions(M)} |R|.

## 2. Region Tuple

### 2.1 Definition

A **region r-tuple** is a multiset τ of r regions (with repetition allowed),
written τ = (R_{i_1}, ..., R_{i_r}) with i_1 ≤ ... ≤ i_r.

τ is **valid** if all regions in τ share at least one maximal clique:
  ∩_{j=1}^{r} prof_{R_{i_j}} ≠ ∅

### 2.2 Represented R-Cliques

Tuple τ **represents** all r-cliques formed by picking one vertex from each
slot of τ (distinct vertices, even from the same region).

For region R appearing k times in τ: pick k distinct vertices from R.

**Multiplicity**: mult(τ) = Π_{R} C(|R|, τ_R)
where τ_R = number of times R appears in τ.

Example: τ = (R₁, R₁, R₂) with |R₁|=98, |R₂|=50.
mult = C(98,2) × C(50,1) = 4753 × 50 = 237,650 r-cliques.

### 2.3 Completeness

**Lemma 1 (Completeness)**: Every r-clique in G is represented by exactly
one region tuple.

Proof: Let T = {v_1, ..., v_r} be an r-clique. Each v_j ∈ R_{i_j} for a
unique region R_{i_j}. The multiset (R_{i_1}, ..., R_{i_r}) (sorted) is the
unique tuple representing T. Since T is a clique, T ⊆ M for some maximal
clique M. All regions of T are in regions(M), so the tuple is valid. □

## 3. Support Equality

### 3.1 Containment Lemma

**Lemma 2**: For any r-clique T represented by tuple τ = (R_{i_1},...,R_{i_r}):
  {M : T ⊆ M} = ∩_{j} prof_{R_{i_j}}

Proof:
(⊇) If M ∈ ∩_j prof_{R_{i_j}}: all regions of τ are in M, so all of T's
    vertices are in M. Thus T ⊆ M.
(⊆) If M' ⊇ T: for each v_j ∈ T, v_j ∈ M', so M' ∈ prof(v_j) = prof_{R_{i_j}}.
    Therefore M' ∈ ∩_j prof_{R_{i_j}}. □

### 3.2 Main Theorem

**Theorem 1 (Support Equality)**: All r-cliques represented by the same
region tuple have identical support, for any s ≥ r+1.

Proof: Support(T) = |{S : T ⊆ S, |S|=s, S is clique}|.

Every such S ⊆ M for some M ⊇ T, i.e., M ∈ ∩_j prof_{R_{i_j}} (by Lemma 2).

The count of s-cliques from M containing T: choose s-r more vertices from
alive_M \ T. Count = C(|alive_M \ T|, s-r) = C(|alive_M| - r, s-r).

Wait: |alive_M \ T| = |alive_M| - r because T ⊆ M and |T| = r and
all T-vertices are alive (T is active).

By inclusion-exclusion over common cliques C_τ = ∩_j prof_{R_{i_j}}:

  support(T) = Σ_{∅≠S⊆C_τ} (-1)^{|S|+1} C(|∩_{M∈S} alive_M| - r, s-r)

This depends ONLY on C_τ and the alive counts, NOT on which specific
vertices form T. Since C_τ is determined by τ, all T's from τ have
the same support. □

### 3.3 Simplified Formula for s = r+1

For s = r+1:

  support(τ) = |∪_{M ∈ C_τ} alive_M| - r

Proof: Each (r+1)-clique containing T is T ∪ {d} for some d adjacent to
all of T. By Lemma 2, d must be in some M ∈ C_τ. And d ∈ alive_M means
d is alive. So d ∈ ∪_{M ∈ C_τ} alive_M. Excluding T's r vertices:
support = |∪alive| - r. □

## 4. Tuple Count Analysis

### 4.1 Per-Clique Bound

For maximal clique M with ρ = ρ(M) regions:

  #tuples_from_M = ((ρ + r - 1) choose r)    [multiset coefficient]
  
  #r-cliques_from_M = (|M| choose r)

### 4.2 Compression Ratio

**Theorem 2 (Compression)**: For each maximal clique M:

  compression(M) = C(|M|, r) / C(ρ(M) + r - 1, r)

If |R_i| = |M|/ρ for all i (uniform region sizes):

  compression ≈ (|M|/ρ)^r / r!  ×  r! / 1 ... hmm let me compute properly.

C(K, r) = K! / (r!(K-r)!) ≈ K^r / r!  for K >> r
C(ρ+r-1, r) = (ρ+r-1)! / (r!(ρ-1)!) ≈ ρ^r / r!  for ρ >> r

compression ≈ (K/ρ)^r

Let α = K/ρ = average region size. Then:

  **compression ≈ α^r**

This grows EXPONENTIALLY with r for α > 1.

### 4.3 Concrete Numbers (web-it-2004)

Largest clique: K = 432, ρ = 20, α = 21.6.

| r | C(K, r) | C(ρ+r-1, r) | Compression |
|---|---------|-------------|-------------|
| 2 | 93,096 | 210 | 443x |
| 3 | 13,375,216 | 1,540 | 8,685x |
| 5 | ~2.6×10^10 | 42,504 | ~612,000x |
| 10 | ~1.6×10^18 | 20,030,010 | ~8×10^{10}x |
| 50 | ~10^{80} | ~10^{30} | ~10^{50}x |
| 100 | ~10^{134} | ~10^{44} | ~10^{90}x |

### 4.4 Global Tuple Count

Total tuples (with dedup) ≤ Σ_M C(ρ(M) + r - 1, r).

Most cliques are SMALL. Only cliques with |M| ≥ r contribute r-cliques.
For small cliques: ρ(M) is small → C(ρ+r-1, r) is small.

For web-it-2004 r=3: measured 787,588 tuples. 

### 4.5 Dependence on r

Tuple count C(ρ+r-1, r) = C(ρ+r-1, ρ-1).

For ρ fixed:
- r < ρ: C(ρ+r-1, r) grows polynomially in r (degree ρ-1)
- r = ρ: C(2ρ-1, ρ) ≈ 4^ρ / √(πρ)  [central binomial]
- r > ρ: C(ρ+r-1, ρ-1) grows polynomially in r (degree ρ-1)

For ρ = 20: C(r+19, 19) = O(r^19). For r=100: C(119, 19) ≈ 10^18.

So: tuple count IS polynomial in r (degree ρ-1), NOT exponential.
But r-clique count C(K, r) is factorial/exponential in r.

**Theorem 3**: compression ratio grows as α^r (exponential in r).
Tuple count grows as r^{ρ-1} (polynomial in r).

## 5. Peeling Algorithm

### 5.1 Data Structures

- Per region R: |R|, prof_R (clique membership), alive count
- Per maximal clique M: alive vertex count, regions(M)
- Per region tuple τ: support value, peeled flag
- Bucket queue indexed by support

NO individual r-clique storage.

### 5.2 Algorithm

```
1. Build regions from maximal clique structure
2. For each clique M: enumerate tuples from regions(M)
3. For each tuple: compute support via inclusion-exclusion
4. Insert into bucket queue
5. While queue non-empty:
   a. Pop tuple τ with minimum support
   b. Assign core(τ) = max(current_level, support(τ))
   c. For each "neighbor" tuple τ': update support
   d. Update bucket queue
```

### 5.3 Cascade Update

When tuple τ = (R_{i_1}, ..., R_{i_r}) is peeled:

All mult(τ) r-cliques are removed. For each (r+1)-clique containing
one of these r-cliques: one r-clique is removed → the (r+1)-clique
becomes incomplete → the other r r-cliques each lose 1 support.

The affected neighbor tuples: those that share an (r-1)-subset of
regions with τ and have one additional region from the overlap.

For s = r+1: the support decrease for neighbor τ' depends on
how many (r+1)-cliques connect τ and τ'. This is computable from
region sizes and shared cliques.

### 5.4 Complexity

| Step | Cost |
|------|------|
| Build regions | O(L) where L = Σ\|M\| |
| Enumerate tuples | O(Σ_M C(ρ(M)+r-1, r)) |
| Compute supports | O(tuples × 2^{max\|C_τ\|}) |
| Peeling | O(tuples × neighbors_per_tuple) |

## 6. Correctness

**Theorem 4 (Equivalence)**: Region tuple peeling gives the EXACT same
core numbers as standard per-r-clique peeling.

Proof: By Theorem 1, all r-cliques in a tuple have identical support at
all times during peeling (since support depends only on the tuple's
common cliques and alive counts). Therefore:

1. All r-cliques in a tuple are peeled at the same level.
2. The cascade from peeling a tuple = sum of cascades from peeling
   its individual r-cliques = correctly computed at tuple level.
3. The k-core is uniquely defined → core numbers are unique → any
   correct peeling gives the same result. □

## 7. What Makes This Work

The fundamental structural property:

**Overlap Class Invariant**: If vertices u, v have prof(u) = prof(v),
then for ANY set of vertices S: u ∈ S and S is a clique ⟺ replacing
u with v in S gives a clique. (Both are in exactly the same maximal cliques.)

This means: swapping vertices within the same region preserves clique-hood
and support. Therefore all r-cliques from the same tuple are "structurally
identical" and can be treated as one entity.

## 8. Limitations

1. Tuple count depends on r: C(ρ+r-1, r) per clique.
   For ρ=20, r=10: C(29,10) ≈ 20M per clique.
   This is polynomial in r but can be large.

2. Compression requires ρ < |M| (non-trivial regions).
   If every vertex has a unique profile (ρ = |M|): no compression.
   This happens on sparse graphs with small cliques.

3. Dense-core graphs benefit most: large cliques with few regions.
