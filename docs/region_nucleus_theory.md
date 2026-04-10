# Region-Based (r,s)-Nucleus Decomposition: Unified Theory

## 1. Definitions

### 1.1 Standard (r,s)-Nucleus Decomposition

- **r-clique**: complete subgraph on r vertices
- **s-clique**: complete subgraph on s vertices (s > r)
- **Support**: deg_s(T) = |{S : T ⊆ S, S is s-clique in current alive graph}|
- **k-(r,s)-core**: maximal subgraph where every r-clique has support ≥ k
- **Core number**: κ(T) = max k such that T is in the k-(r,s)-core
- Core values are **non-decreasing** in peeling order
- The k-core is **uniquely defined** for each k

### 1.2 Maximal Cliques and Regions

- **Maximal clique set**: M = {M_1, ..., M_m}
- **Profile**: prof(v) = {M ∈ M : v ∈ M}
- **Region**: R = maximal set of vertices with identical profile
  - Region set: R = {R_1, ..., R_p}, p ≤ n
  - Every vertex belongs to exactly one region

## 2. Proven Lemmas

### Lemma 1 (Region Clique Property)
All vertices of a region R form a clique.

*Proof*: Any u, v ∈ R share at least one maximal clique M ∈ prof_R. □

### Lemma 2 (Containment)
For any r-clique T ⊆ R: {M : T ⊆ M} = prof_R.

*Proof*:
- (⊇) T ⊆ R ⊆ M for all M ∈ prof_R.
- (⊆) If M' ⊇ T, pick t ∈ T ⊆ R. Then M' ∈ prof(t) = prof_R. □

### Theorem 1 (Same-Region Support Equality)
All same-region r-cliques from R have identical support, for any s.

*Proof*: By Lemma 2, every T ⊆ R has containing maximal cliques = prof_R.
The s-cliques containing T depend only on prof_R and alive counts. □

### Theorem 2 (Cross-Region Support Bound)
For r-clique T with vertices from regions R_1, ..., R_k:
support(T) ≤ σ(R_i) for each i.

*Proof*: {M : T ⊆ M} = ∩_i prof_{R_i} ⊆ prof_{R_i}. □

## 3. Region-Only Peeling: UPPER BOUND (Not Exact)

### 3.1 Definition

σ(R) = support of any same-region r-clique from R.

For s = r + 1:
  σ(R) = |∪{M ∈ prof_R} alive_M| - r

For general s:
  σ(R) = Σ_{∅≠S⊆prof_R} (-1)^{|S|+1} C(|∩_{M∈S} alive_M| - r, s - r)

### 3.2 Region Peeling is NOT Exact for r ≥ 2

**Counterexample** (r=2, s=3, 7 vertices):

Maximal cliques: {0,1,2,3}, {0,1,4}, {2,3,5}, {4,5,6}.
Region R_A = {0,1}: prof = {{0,1,2,3}, {0,1,4}}.

- σ(R_A) = |{0,1,2,3,4}| - 2 = 3
- Cross-region edge {0,4}: support = |{0,1,4}| - 2 = 1
- Per-edge peeling: {0,4} peeled at level 1, cascade reduces σ(R_A) to 2
- Region peeling: R_A peeled at level 3 (cross-region cascade delayed)

**Result**: region_core(R_A) = 3, true_core({0,1}) = 2. Overestimate.

### 3.3 Why Region Peeling Overestimates

Region peeling does NOT account for the early peeling of cross-region r-cliques.
These have support < σ(R), and their cascade reduces σ before R is peeled.
Region peeling delays this cascade → σ stays inflated → overestimated cores.

### 3.4 Region Peeling IS Exact for r = 1

For r = 1: every r-clique (= vertex) is same-region. No cross-region r-cliques
exist. So the cascade delay issue does not arise. Exact. □

### 3.5 Upper Bound Property

**Theorem 3**: region_core(R) ≥ true_core(T) for all same-region T ⊆ R.

*Proof sketch*: The delayed cascade keeps σ(R) ≥ the true per-r-clique support.
Since we peel at a higher level, the core value is at least as high. □

## 4. The Unified Theory: Region-Tuple Decomposition

### 4.1 Key Insight

The problem with region-only peeling: it ignores cross-region r-cliques.
The fix: track **region tuples** instead of individual r-cliques.

### 4.2 Definition

**Region r-tuple**: an ordered multiset (R_{i_1}, ..., R_{i_r}) with
i_1 ≤ ... ≤ i_r, where all r regions share at least one common maximal clique
(∩_j prof_{R_{i_j}} ≠ ∅).

Each region tuple represents a GROUP of r-cliques:
- Pick one vertex from each region in the tuple
- The number of represented r-cliques depends on the tuple type:
  - All distinct regions: |R_{i_1}| × ... × |R_{i_r}|
  - With repetitions: multinomial coefficient

### 4.3 Fundamental Theorem

**Theorem 4 (Tuple Support Equality)**: All r-cliques represented by the same
region tuple have identical support.

*Proof*: Let tuple τ = (R_{i_1}, ..., R_{i_r}). For any r-clique T represented
by τ: {M : T ⊆ M} = ∩_j prof_{R_{i_j}} (by Lemma 2 applied to each region).
This set depends only on the tuple, not the specific vertices chosen. □

**Corollary**: Region-tuple peeling is EQUIVALENT to per-r-clique peeling.
Since all r-cliques in a tuple have the same support, peeling a tuple =
peeling all its r-cliques simultaneously = peeling them one by one (same level).

### 4.4 Tuple Support Formula

For tuple τ = (R_{i_1}, ..., R_{i_r}) with common cliques C_τ = ∩_j prof_{R_{i_j}}:

**For s = r + 1:**
  support(τ) = |∪{M ∈ C_τ} alive_M| - r

**For general s:**
  support(τ) = Σ_{∅≠S⊆C_τ} (-1)^{|S|+1} C(|∩_{M∈S} alive_M| - r, s - r)

### 4.5 Tuple Count: Always Less Than R-Clique Count

**Theorem 5 (Compression)**: For any maximal clique M, let ρ(M) = number of
regions in M. Then:

  #tuples_from_M = C(ρ(M) + r - 1, r)  (with repetition)
  #r-cliques_from_M = C(|M|, r)

Since ρ(M) ≤ |M| and C(K, r) is increasing in K:
  #tuples_from_M ≤ #r-cliques_from_M

And strict inequality when any region in M has size > 1:
  ρ(M) < |M| ⟹ C(ρ(M)+r-1, r) << C(|M|, r)

**Compression ratio per clique**: C(|M|, r) / C(ρ(M)+r-1, r).
For |M| = 432, ρ(M) = 20, r = 3: C(432,3) / C(22,3) = 13,375,216 / 1,540 ≈ 8,685x
For |M| = 432, ρ(M) = 20, r = 10: C(432,10) / C(29,10) ≈ 10^{18} / 10^7 ≈ 10^{11}x

### 4.6 Cascade Update for Tuples

When tuple τ_peeled is peeled, for each affected tuple τ':

The support of τ' decreases by the number of COMPLETE s-cliques that
(a) contained an r-clique from τ_peeled, AND
(b) contain an r-clique from τ', AND
(c) all OTHER r-cliques of the s-clique are still alive (unpeeled).

For s = r + 1: each s-clique has r+1 r-cliques. We need:
- One r-clique from τ_peeled (just peeled)
- One r-clique from τ' (the affected tuple)
- All r-1 other r-cliques of the s-clique are in UNPEELED tuples

This is computable at the tuple level:
- The s-clique's r+1 r-cliques correspond to r+1 tuples
- Check that the other r-1 tuples (besides τ_peeled and τ') are unpeeled

For general s: each s-clique has C(s,r) r-cliques. The cascade involves
checking that all C(s,r) - 2 other r-cliques (tuples) are unpeeled.

### 4.7 Data Structures

**Stored per tuple**: tuple key (r region IDs), support value, peeled flag.
**Stored per region**: profile, alive count.
**Stored per maximal clique**: alive count, list of regions.
**NO individual r-clique data**.

### 4.8 Complexity

| Step | Cost |
|------|------|
| Build regions | O(L) where L = Σ\|M_i\| |
| Enumerate tuples | O(Σ_M C(ρ(M)+r-1, r)) |
| Compute support | O(tuples × 2^{max\|prof\|}) |
| Cascade peeling | O(tuples × neighbors × r) |

**Key**: #tuples << #r-cliques on dense graphs. The improvement ratio
is C(|M|,r) / C(ρ(M)+r-1,r) per clique, which is MASSIVE when ρ(M) << |M|.

## 5. Why This Is a Unified Theory

### 5.1 Works for ALL r and s

- The tuple support formula (Section 4.4) applies for any r, s.
- The cascade (Section 4.6) applies for any r, s.
- The compression bound (Section 4.5) holds for any r.

### 5.2 Always At Least As Fast As Standard

#tuples ≤ #r-cliques always (Theorem 5). Equal only when every region
has size 1 (worst case: every vertex has a unique maximal clique profile).

### 5.3 Improvement Grows with Clique Size

For a maximal clique M with ρ(M) regions:
  compression = C(|M|, r) / C(ρ(M)+r-1, r)

This grows SUPER-EXPONENTIALLY with |M| for fixed ρ(M) and r.

### 5.4 Exactness

Region-tuple peeling is EXACT: it gives identical core numbers to
per-r-clique peeling (Theorem 4). No approximation.

## 6. Experimental Validation

| Graph | r | s | r-cliques | Tuples | Compression | Time |
|-------|---|---|-----------|--------|-------------|------|
| dblp-core30 | 3 | 4 | 686K | 1,706 | 402x | 2ms |
| web-it-2004 | 3 | 4 | 338M | 787K | 430x | 1.75s |
| dblp-core30 max core | 3 | 4 | 111 (ST) | 111 (tuple) | exact | |
| web-it-2004 max core | 3 | 4 | timeout | 429 | | |

## 7. Open Questions

### 7.1 Can We Achieve r-Independent Complexity?

Region-only peeling is O(L·R_max), independent of r, but OVERESTIMATES for r≥2.
Region-tuple peeling is EXACT but has O(Σ C(ρ(M)+r-1, r)) tuples.

Question: is there a way to get EXACT results with r-independent complexity?

Current answer: seems unlikely. Cross-region r-cliques inherently create
r-dependent interactions. But the tuple count is much smaller than the
r-clique count, making the approach practical.

### 7.2 Cascade for General s

For s > r + 1: the cascade involves checking C(s,r) - 2 other tuples per
s-clique, which is more expensive. Need efficient implementation.

### 7.3 Maximal Clique Tagging

Current SDCT_MaxClique tags ~80-100% of maximal cliques correctly.
Incorrect tagging causes ~0.03% of triangles to be untracked.
Need to fix the X propagation in BK for complete coverage.
