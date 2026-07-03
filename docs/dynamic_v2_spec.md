# Dynamic (1,s)-core Maintenance — v2 Specification

**Status**: authoritative design + theory document for the v2 rewrite.
**Audience**: an implementing agent with access to this repository. Read this
document END TO END before writing code. Every design decision here is
load-bearing; several were forced by real bugs found during v0/v1 (see §16).
**Reference implementation of v1**: `src/dynamic_1s_core.cpp` (commit
`131bf62`), harness `bench_dynamic_insert.py`, locality data
`dynamic_locality_out/`.

---

## 1. Problem statement and notation

Undirected simple graph `G = (V, E)`, `n = |V|`. Fix integer `s >= 3`
(s = 2 degenerates to k-core; not our target but everything below still
holds). An **s-clique** is a set of `s` pairwise-adjacent vertices.

For `x ∈ T ⊆ V`, define the **support**
`sup_T(x) = #{ s-cliques K in the induced subgraph G[T] with x ∈ K }`.

Note: `T ⊆ T'  ⟹  sup_T(x) <= sup_T'(x)` (monotone in `T`).

**ℓ-core**: `C_ℓ(G)` = the maximum `T ⊆ V` with `sup_T(x) >= ℓ` for every
`x ∈ T`. (Maximum exists and is unique: if `T1, T2` both satisfy the
property, so does `T1 ∪ T2`, by support monotonicity. Also
`C_{ℓ+1}(G) ⊆ C_ℓ(G)`.) `C_0(G) = V`.

**Core number**: `core_G(x) = max{ ℓ : x ∈ C_ℓ(G) }`.

**Update model**: a single edge insertion `G' = G + e`, `e = (u,v)`,
`e ∉ E`. Given `core_G[·]` (array) and adjacency of `G`, produce
`core_G'[·]` exactly. (Deletion: §13. Batches: §17/M5.)

Shorthand: `c(x) = core_G(x)`, `c'(x) = core_G'(x)`,
`W = N_G(u) ∩ N_G(v)` (common neighborhood, identical in `G'`),
`R* = { x : c'(x) ≠ c(x) }` (the true changed set).

Supports and cores are **integers** but can be astronomically large
(binomially many cliques); store them in `double` (exact for values
< 2^53, same convention as the whole codebase; all counting must use the
same arithmetic as the reference so comparisons are bit-stable).

---

## 2. Peel semantics = core numbers (Lemma 0)

The reference algorithm (V3 and our scoped peel) uses this loop:

```
alive := V; m := 0
while alive nonempty:
    x := argmin_{y ∈ alive} sup_alive(y)     // any tie order
    m := max(m, sup_alive(x))                 // the "max-rule"
    core(x) := m
    alive := alive \ {x}
```

**Lemma 0.** The peel assigns `core(x) = core_G(x)` for every `x`,
regardless of tie order.

*Proof.* Let `A_ℓ = {x : peel value >= ℓ}`. Because `m` is nondecreasing,
there is a first pop `t*` at which `m >= ℓ`; exactly the pops from `t*`
onward get value `>= ℓ`, so `A_ℓ` = the alive set just before `t*`. At
`t*`, `m_old < ℓ` and `m_new = max(m_old, sup(x_{t*})) >= ℓ` forces
`sup_alive(x_{t*}) >= ℓ`; since `x_{t*}` is the minimum, every `y ∈ A_ℓ`
has `sup_{A_ℓ}(y) >= ℓ`. Hence `A_ℓ` is ℓ-satisfying, so `A_ℓ ⊆ C_ℓ(G)`,
giving peel <= core. Conversely, while all of `C_ℓ(G)` is alive, any
`y ∈ C_ℓ(G)` has `sup_alive(y) >= sup_{C_ℓ}(y) >= ℓ`; so the first pop of
a `C_ℓ(G)` member has key `>= ℓ`, which pushes `m >= ℓ` at that pop —
i.e. no `C_ℓ(G)` member is ever assigned a value < ℓ. Hence core <= peel.
∎

Two consequences used throughout:
- **Tie order never affects final values.** (Batching equal-level
  removals is safe.)
- The max-rule matters: pops after `m` rises can have keys `< m` and
  still receive value `m`.

---

## 3. Insertion monotonicity (Lemma 0')

**Lemma 0'.** `G ⊆ G' ⟹ core_G(x) <= core_G'(x)` for all `x`.
*Proof.* `C_ℓ(G)` is ℓ-satisfying inside `G'` too (supports only gain),
so `C_ℓ(G) ⊆ C_ℓ(G')`. ∎

So on insertion `R* = { x : c'(x) > c(x) }` (**risers only**). Any vertex
whose core "drops" in an implementation is a bug (keep the v1 WARN-drop
guard).

---

## 4. Scoped exactness — the Phase-3 foundation

**Definition (conditioned core).** Given a region `R ⊆ V` and boundary
values `b : V\R → ℕ`, for `ℓ >= 1` let `D_ℓ(R,b)` be the maximum
`D ⊆ R` such that every `x ∈ D` has at least `ℓ` s-cliques of `G'`
containing `x` inside `D ∪ { w ∉ R : b(w) >= ℓ }`. (Maximum exists by
union-closure as in §1; `D_ℓ` is nested decreasing in `ℓ` because the
boundary set `{b >= ℓ}` shrinks.) Define
`ccore_{R,b}(x) = max{ ℓ : x ∈ D_ℓ(R,b) }`.

**Pinned peel on (R, b)** (this is what v1 implements): universe
`U = R ∪ ∂R`, `∂R = N(R) \ R`. Region vertices carry live keys
`sup` (within the current alive set, in `G'`); pinned vertices `w ∈ ∂R`
carry the FIXED key `b(w)`. Pop the global minimum key (pinned pops are
just removals; region pops receive `core := m` under the max-rule).
Note every s-clique through `x ∈ R` lies inside `{x} ∪ N(x) ⊆ U`, so `U`
suffices for all region support counting.

**Lemma 1 (pinned peel computes ccore).** For every `x ∈ R` the pinned
peel assigns `ccore_{R,b}(x)`.

*Proof.* Same two directions as Lemma 0 with one extra observation about
the alive set at the first pop with `m >= ℓ` ("the crossing"):
(a) all pinned `w` with `b(w) < ℓ` are already removed at the crossing —
otherwise the min-rule would have selected `w` (its key `b(w) < ℓ <=`
crossing key) before the crossing pop; (b) all pinned with `b >= ℓ` are
still alive — popping one sets `m >= ℓ`, contradicting `m_old < ℓ`.
So at the crossing, alive `= {x ∈ R : final >= ℓ} ∪ {w : b(w) >= ℓ}`,
and the crossing key `>= ℓ` bounds all alive region supports from below.
Direction 1: that alive region set is inside `D_ℓ(R,b)`. Direction 2:
while `m < ℓ`, all of `D_ℓ(R,b)` and all pinned `{b >= ℓ}` are alive by
induction, so any `x ∈ D_ℓ` has `sup_alive(x) >= ℓ`; hence a `D_ℓ`
member popped while `m < ℓ` would have key `>= ℓ` and push `m >= ℓ` at
its own pop — i.e. it still receives a value `>= ℓ`. ∎

**Lemma 2 (scoped exactness).** If `R ⊇ R*` (i.e. `c'(y) = c(y)` for all
`y ∉ R`) and `b(w) = c(w)` on `∂R`, then `ccore_{R,b}(x) = c'(x)` for
every `x ∈ R`.

*Proof.* (>=) `D*_ℓ := C_ℓ(G') ∩ R` is a valid candidate for `D_ℓ(R,b)`:
a witness member `z ∉ R` has `c'(z) >= ℓ` and `c'(z) = c(z) = b(z)`, so
it lies in the boundary set. Hence `D*_ℓ ⊆ D_ℓ(R,b)`.
(<=) `T_ℓ := D_ℓ(R,b) ∪ (C_ℓ(G') \ R)` is ℓ-satisfying in `G'`:
for `z ∈ D_ℓ` its conditioned witnesses live in
`D_ℓ ∪ {w ∉ R : c(w) >= ℓ}`; `c(w) = c'(w) >= ℓ` puts those boundary
vertices in `C_ℓ(G') \ R ⊆ T_ℓ`. For `z ∈ C_ℓ(G')\R` its witnesses live
in `C_ℓ(G') ⊆ T_ℓ` (using `C_ℓ(G') ∩ R = D*_ℓ ⊆ D_ℓ`). Hence
`T_ℓ ⊆ C_ℓ(G')`. ∎

**Consequence — the whole correctness burden reduces to one job:
produce any `R ⊇ R*`.** Then ONE pinned peel finishes the update. All of
v1's expansion rounds, riser windows, and promotion gates existed only
because v1 could not certify `R ⊇ R*` up front. v2 fixes exactly this.

---

## 5. Structure of the change set

Write `D_ℓ = C_ℓ(G') \ C_ℓ(G)` (new survivors at level ℓ). Facts:
`x ∈ D_ℓ ⟺ c(x) < ℓ <= c'(x)`; `R* = ∪_ℓ D_ℓ`;
`C_ℓ(G') \ D_ℓ = C_ℓ(G)`.

A **witness of `z` at level ℓ** is an s-clique `K ∋ z` of `G'` with
`K ⊆ C_ℓ(G')`. Every `z ∈ C_ℓ(G')` has `>= ℓ` witnesses at level ℓ.

**Lemma 3 (connectivity of new survivors — set version).**
Let `Z ⊆ D_ℓ` be nonempty and suppose that every witness (at level ℓ) of
every `z ∈ Z` (i) does not contain both `u` and `v`, and (ii) satisfies
`K ∩ D_ℓ ⊆ Z`. Then contradiction — i.e. every such `Z` must contain a
member with a witness containing `{u,v}`, or a witness touching
`D_ℓ \ Z`.

*Proof.* Under (i)+(ii), each witness `K` of `z ∈ Z` avoids `e`, hence
is a clique of `G`, and `K ⊆ (C_ℓ(G') \ D_ℓ) ∪ Z = C_ℓ(G) ∪ Z`. Then
`T' = C_ℓ(G) ∪ Z` is ℓ-satisfying **in G** (old members are supported
within `C_ℓ(G)` alone; `Z` members by the above). So `Z ⊆ C_ℓ(G)`,
contradicting `Z ⊆ D_ℓ`. ∎

**Corollary 3a (chains).** Define the witness graph on `D_ℓ`: join
`z ~ z'` when some witness of `z` (at level ℓ) contains `z'`. Every
connected component of this graph contains a vertex having a witness
that contains both `u` and `v`. Hence every `x ∈ D_ℓ` is linked to the
insertion by a finite chain `x = x_0 ~ x_1 ~ … ~ x_m`, all `x_i ∈ D_ℓ`,
where `x_m` has a witness `K ∋ u, v`. Chain-adjacent vertices share a
clique, hence are **graph-adjacent**. Members of a `{u,v}`-witness other
than `u,v` lie in `W`.

**Corollary 3b (levels of chain members).** All chain members satisfy
`c(x_i) < ℓ <= c'(x_i)`, so each `x_i ∈ R*`; and since chains exist at
EVERY level `ℓ ∈ (c(x), c'(x)]`, in particular at `ℓ_x := c(x)+1`.

**Corollary 3c (no new cliques ⟹ no change).** If `G[W]` contains no
(s−2)-clique, then no s-clique of `G'` contains both `u,v`, every
`C_ℓ(G')`-witness is a `G`-clique, and Lemma 3 with `Z = D_ℓ` forces
`D_ℓ = ∅` for all ℓ. **Early exit** (this covers the 16–55% of
zero-change inserts measured in `dynamic_locality_out/`).

---

## 6. Phase 2 — candidate discovery (the heart of v2)

### 6.1 The static admission test

For any vertex `z` let `TS(z) = sup^{G'}_V(z)` = total s-clique support
of `z` in `G'` (post-insertion adjacency, no filter). Note
`c'(z) <= TS(z)` (support within any set is at most total support).

**Definition (admission test).** For a vertex `y`, with `ℓ_y = c(y)+1`:

```
OS(y) = #{ s-cliques K ∋ y in G' : every z ∈ K\{y} has
           c(z) >= ℓ_y   OR   TS(z) >= ℓ_y }
PASS(y)  ⟺  OS(y) >= ℓ_y
```

Because `c(z) <= TS(z)` always, the first disjunct is implied by the
second; operationally you evaluate `c(z) >= ℓ_y` first (free) and
compute `TS(z)` lazily (one SCT over `N(z)`, memoized per update) only
for members failing it.

**Crucial property: the test is STATIC — it does not depend on the
current candidate set.** No retesting, no worklist fixpoint; each vertex
is tested at most once per update.

**Lemma 4 (test completeness).** Every `y ∈ R*` passes:
its `>= ℓ_y` witnesses at level `ℓ_y` (which exist since
`c'(y) >= ℓ_y`) have members `z` with `c'(z) >= ℓ_y`, hence
`TS(z) >= c'(z) >= ℓ_y`. ∎

*(Why the test cannot be strengthened naively: mutual-rise groups are
real — v1's "zero-net-rise deadlock" was a 6-vertex ring at core 55
whose members support each other; any admission rule that requires a
full witness with already-admitted/high members deadlocks on rings.
The TS-disjunct is exactly what admits ring members without an order.
Rings can also span MIXED levels (a c=53 vertex and a c=55 vertex rising
together to 56 is legal by Corollary 3b), so same-level-only optimism —
what v1's promotion gate used — is provably insufficient in general,
even though it happened to pass our test suite.)*

### 6.2 The trigger closure

```
C := { u, v }                     // unconditional seeds
frontier := { u, v }
tested := { u, v }
while frontier nonempty:
    pop x from frontier
    for y ∈ N_G'(x), y ∉ tested:
        tested += y
        if PASS(y):  C += y; frontier += y
```

(Equivalently: `C` = the connected component of `{u,v}` in the subgraph
induced by PASS-ing vertices, unioned with `{u,v}`.)

**Theorem 5 (discovery completeness).** `R* ⊆ C`.

*Proof.* Let `y ∈ R*` and take its chain at level `ℓ_y = c(y)+1`
(Corollary 3a/3b): `y = x_0 ~ … ~ x_m`, `x_m` sharing a witness with
`u, v`. Induct backwards from the seed: `x_m` is graph-adjacent to `u`
(co-members of one clique), so it gets tested once `u ∈ C` (immediately);
`x_m ∈ R*` passes by Lemma 4, so `x_m ∈ C` and enters the frontier.
Then `x_{m−1}` is graph-adjacent to `x_m`, gets tested, passes
(∈ `R*`), … down to `x_0 = y`. Each step only needs the previous chain
member to be IN `C` (for the trigger) plus the static test (no
dependence on how much of the ring is admitted yet). ∎

**Remark (why not test only clique-sharing neighbors):** the trigger
uses plain graph adjacency `N(x)`, a superset of clique-sharing. This
only adds tests (safe); restricting the trigger to "shares an s-clique
with x" is also sound (chain partners share cliques) and cheaper per
frontier vertex when implemented as: test `y` only if
`|N(y) ∩ N(x)| >= s−2`. Optional micro-optimization; keep plain `N(x)`
in the first correct version.

### 6.3 Eviction (tightening, C-dependent)

After the closure completes (never interleaved — see below), prune:

```
EOS_C(x) = #{ s-cliques K ∋ x in G' : every z ∈ K\{x} has
              c(z) >= c(x)+1   OR   (z ∈ C AND TS(z) >= c(x)+1) }
repeat until stable:
    if some x ∈ C \ {u,v} has EOS_C(x) < c(x)+1:  C := C \ {x}
```

(`u,v` stay: they are the insertion endpoints and must be in the peel
region; they cost nothing extra.)

**Lemma 6 (eviction safety).** If `R* ⊆ C` before eviction, no member of
`R*` is ever evicted (so `R* ⊆ C` is an invariant and the final `C`
still satisfies Lemma 2's hypothesis).

*Proof.* Induction over evictions. For `x ∈ R*`, its level-`ℓ_x`
witnesses have members `z` with `c'(z) >= ℓ_x`; either `c(z) >= ℓ_x`
(first disjunct) or `z ∈ R* ⊆ C` (induction) with
`TS(z) >= c'(z) >= ℓ_x` (second disjunct). So `EOS_C(x) >= ℓ_x`, and
`x` is never selected. ∎

When `x` is evicted, only neighbors of `x` can lose EOS value; update
their `EOS` by delta subtraction (count cliques through the pair
`(neighbor, x)` under the same member-filter — the `held = 2` SCT of
§11) or recompute; cascade until stable.

**HARD CONSTRAINT — strict phase order.** Eviction must start only
after the closure has fully saturated. Interleaving is UNSOUND: Lemma
6's proof needs `R* ⊆ C` at eviction time; evicting early, while a ring
member's partners are not yet admitted, can kill a true riser. An
implementing agent will be tempted to interleave for speed. DO NOT.

### 6.4 What Phase 2 replaces

v2 has **no expansion rounds, no riser window, no same-level promotion
gate, no 40-round cap**. Phase 2's output `C` (plus Theorem 5 + Lemma 6
+ Lemma 2) certifies that ONE pinned peel suffices. This removes v1's
dominant cost (up to 38 full re-peels; profiled rounds paying 850 ms to
add one vertex).

### 6.5 Size control and principled fallback

`C` can exceed `R*` (false positives are pruned partly by eviction, and
harmlessly resolved by the peel — a non-riser candidate just re-derives
its old core). On adversarial shells (dblp s=3: v1 saw candidate regions
up to 13,755) `C` may still be large; each closure step costs one static
test (local SCT), not a re-peel, so even the crawl case is ~10^4 cheap
tests. Nevertheless implement a **cap**: if `|C| > CAP` (default:
`max(4096, n/16)`) abort the incremental path, print `FALLBACK`, and let
the driver run a full static recompute. This is the honest hybrid; count
fallbacks in stats. (Tightening knob if measurements demand it: replace
`TS(z)` in the tests by any cheaper-but-still-`>= c'(z)` upper bound, or
add leveled candidate caps `UB(z) = OS(z)` — sound since
`c'(z) <= OS(z)` by Lemma 4's argument applied at every level up to
`c'(z)`. Only pursue after measuring.)

---

## 7. Phase 1 — seed computation and early exit

1. Compute `W = N(u) ∩ N(v)` (sorted-merge).
2. Count (s−2)-cliques in `G[W]` with ONE local SCT (closed forms, §11).
   If the count is zero: **no core changes** (Corollary 3c) — output no
   CHANGED lines, print STATS, exit. This must be the fast path
   (microseconds): do it before any allocation proportional to `n`
   beyond what graph loading already did.
3. Insert `e` into the adjacency NOW (before Phase 2 — all Phase-2/3
   counting is in `G'`).

Per-vertex new-clique counts `Δ(x)` are NOT needed for correctness in
v2 (the static test on `W`-members sees the new cliques automatically,
because tests count in `G'`). Emit the total as a stat only.

---

## 8. Phase 3 — one pinned peel

Run the v1 pinned peel ONCE on region `C` (final, post-eviction):

- Universe `U = C ∪ ∂C`, `∂C = N(C)\C`; pinned keys `b(w) = c(w)`.
- Pinned sorted ascending by `c`, consumed by cursor (never in the
  scan pool). Region pool: linear scan-min is fine (|C| small); switch
  to a bucket structure only if profiling demands.
- **Max-rule** `m := max(m, key)` exactly as §2; region pops get
  `newCore := m`.
- All three v1 exact optimizations carry over verbatim, with proofs:

**L0 fast-start.** `L0 = min_{x∈C} c(x)`; initialize `m := L0`; pinned
with `c < L0` are never alive (excluded from `U`). *Proof*: every
region vertex's initial key `>= c(x) >= L0` (its `c(x)`-level witnesses
in `G` have members with `c >= c(x) >= L0`, all alive at init), so no
region pop can occur below `L0`; pinned below `L0` would all be popped
first and only ratchet `m` up to at most `L0`, and their removals do not
change any region key that matters (their cliques are exactly those
excluded by the τ-view below). Identical final values.

**Per-vertex τ-view.** Region vertex `x` accounts ONLY cliques whose
pinned members all have `c >= τ_x := c(x)` (region members always
count): in `supportOf(x)` and in `supportDelta(x, z)` filter members by
`inView(y, τ_x) ⟺ alive[y] AND (inRegion[y] OR c(y) >= τ_x)`; skip the
delta call entirely for pinned `z` with `c(z) < τ_x`.
*Proof of exactness*: (i) the τ-view initial key still satisfies
`key(x) >= c(x)` (same witness argument as L0 — those witnesses' pinned
members have `c >= c(x) = τ_x`, so the view keeps them); therefore `x`
is never the scan-minimum while `m < τ_x`: any pinned with key
`< τ_x <= key(x)` is popped first by the min-rule. (ii) By the time
`m >= τ_x`, every pinned below `τ_x` has been drained, so from the first
moment `x` can pop, its τ-view key equals the unfiltered key. Pop
decisions coincide at every step; all outputs identical. *(This removed
the 10k-delta hub storms in v1: pinned below τ generate no work for x.)*

**Exact delta subtraction.** When removing any `z` (pinned or region
pop), for each alive region neighbor `y` of `z` (respecting the τ-skip):
`key[y] -= supportDelta(y, z)` computed BEFORE `alive[z] := 0`; the
delta counts s-cliques containing both `y, z` with all other members in
`y`'s view (`held = 2` SCT on `N(y) ∩ N(z)` filtered, §11). Keys stay
exact at every step. **Never batch removals across distinct levels with
deferred recounts** — see §16 scar (b).

Guards (keep, cheap): negative key ⟹ WARN + clamp 0; final
`newCore < c` ⟹ WARN drop (both indicate bugs; they fired during v1
development and caught real errors).

Output: `CHANGED x old new` for every `x ∈ C` with `newCore ≠ c(x)`,
then a `STATS` line (§14).

---

## 9. Full update algorithm (assembled)

```
insert(u, v):
  Phase 1: W; count new cliques; if 0 → STATS; exit.   // §7
           insert e into adjacency (G' from here on).
  Phase 2: closure from {u,v} with static PASS test     // §6.2
           (memoize TS; tested-set prevents re-tests);
           if |C| > CAP → FALLBACK; exit.
           eviction cascade to fixpoint.                // §6.3 (strictly after)
  Phase 3: one pinned peel on C (τ-view, L0, deltas).   // §8
  Emit CHANGED + STATS.
```

Persistent state between updates: `core[]` and adjacency ONLY. This is a
selling point (no index maintenance); do not add hidden persistent
caches without flagging it (support caches are a legitimate future
optimization but change the memory story).

---

## 10. Complexity (honest statement)

Let `T = tested set ⊆ N[C]`, `T_cnt(y)` = cost of one filtered SCT on
`N(y)`-scale sets. Per update:
`O( Σ_{y∈T} T_cnt(y)  +  eviction cascades  +  PinnedPeel(C) )`,
fully instance-sensitive; with the cap, worst case = one static
recompute (never asymptotically worse than the baseline). Do NOT claim
`|C| = O(|R*|)` — false in general (shells). Report measured `|C|/|R*|`
distributions instead.

---

## 11. Counting engine

All phases use one primitive: count s-cliques through given held
vertices with a member predicate. Current engine (v1, correct,
reuse-then-optimize): pivoter recursion `sctCount(P, held, piv)` in
`src/dynamic_1s_core.cpp`, with invariants:

- Node state: held set size `held` (vertices REQUIRED in the clique),
  candidate list `P` (all adjacent to every held and every pivot
  vertex), pivot count `piv` (accumulated Π; all Π vertices pairwise
  adjacent and adjacent to all of P and all held — maintained by the
  recursion, relied upon by the closed forms).
- Leaf `P = ∅`: contributes `C(piv, s − held)`.
- Recursion: choose pivot `p ∈ P` (maximizing `|P ∩ N(p)|`; ANY choice
  is correct, the max only shrinks the tree); branch 1 (pivot):
  `(P ∩ N(p), held, piv+1)`; then remove `p` from the pool and for each
  `v ∈ pool \ N(p)` in order: branch `(pool ∩ N(v), held+1, piv)` then
  remove `v` from the pool. (No clique contains both `p` and a
  non-neighbor of `p`, which is why `p` is deleted before the hold
  branches — do not "fix" this.)
- Exact closed-form base cases (proved by induction in v1, keep both):
  `held == s−1 → piv + |P|`;
  `held == s−2 → C(piv,2) + piv·|P| + E(G[P])` where `E` = #edges
  inside `P`.
- Per-vertex attribution at a leaf `(H, Π)` counting k-cliques (used in
  Phase 1 on `G[W]` with `k = s−2`): a held vertex gets
  `C(|Π|, k−|H|)`; a pivot vertex gets `C(|Π|−1, k−|H|−1)`.
- `nCr` via Pascal's triangle in double (matches repo convention).

**M1 optimization (bit-identical, gate on the harness):** build a local
compact index for the update's working set (`U = C ∪ ∂C` in Phase 3;
`N[y]`-scoped for Phase 2 tests): map to dense local ids, adjacency rows
as bitsets (`ceil(|U|/64)` words); `P` as bitset; pivot selection =
`popcount(P & row(z))`; branch sets = `P & row(v)`; member predicates
(alive / inRegion / c >= τ / TS >= ℓ) folded into precomputed row masks.
Expected 10–100× constant-factor gain over the current per-call
`std::vector` merge-intersects. MUST produce bit-identical outputs (same
counts ⟹ same peel ⟹ same cores).

---

## 12. Verification protocol (non-negotiable)

Two tiers, per the repo's standing rule: correctness before performance,
per-VERTEX diffs (a distribution check once hid a real swap bug — §16).

- **Tier 1 (exact, per edge)**: `bench_dynamic_insert.py` — for each
  sampled edge `e`: reference `core(G−e)` and `core(G)` from the V3
  binary (original-id space via `PIVOTER_DUMP_MAPPING`); run the
  prototype insert on `(G−e) + e`; require merged result ==
  `core(G)` on every vertex (missing id ⟺ core 0). GATE: **0 mismatches
  on all 4 configs × 300 edges** (com-dblp s=3,5; soc-Epinions1 s=3,5;
  seed 42 — the standing 1188-edge gate; 12 dblp-s3 edges may FALLBACK
  by cap, which is acceptable but must be counted and reported).
- **Tier 2 (streaming, big graphs)**: on soc-pokec (s=3,4): remove k=50
  random edges from `G`, compute `core(base)` once, insert them back
  one at a time through the prototype, verify the FINAL state equals
  `core(G)` per-vertex (one reference run total). Any mid-stream error
  surfaces as a final mismatch.

Report failures verbatim; never average away a mismatch.

## 13. Performance acceptance criteria (the bar that matters)

Baseline = **V3 peel-only** (`STV3_PEEL_US`), NOT build+peel — the CPI
build is a one-time cost the dynamic algorithm does not pay. Current
peel-only numbers (this machine, single-thread): com-dblp s=3: 32.6 ms,
s=5: 12.6 ms; soc-Epinions1 s=3: 35.3 ms, s=5: 61.0 ms.

Targets for v2 (measured on the Tier-1 sweeps):
- median insert time: **>= 10× faster than peel-only** on all 4 configs;
- p90: **not slower than peel-only** on any config;
- zero fallbacks on Epinions; fallbacks on dblp-s3 reported and few;
- update cost must scale with `|C|`, not with rounds (there are no
  rounds); attach the `DYN1S_DEBUG`-style phase breakdown
  (test time / eviction time / peel time) for the slowest 5 edges of
  each config.

## 14. Engineering contract

- CLI unchanged: `dynamic_1s_core <base.edges> <s> <core_base.tsv> <u> <v>`;
  output `CHANGED x old new` lines + one `STATS` line; on cap:
  `FALLBACK` line + `STATS`, no CHANGED lines. Extend STATS to:
  `region= pinned= tested= evicted= l0= pinned_skipped= newcliques= insert_us=`
  plus phase timings `p1_us= p2_us= p3_us=`.
- Timing starts at edge insertion into adjacency (exclude file loads).
- Build: add to existing CMake target `dynamic_1s_core`; **never build
  with more than `-j 12`** (repo hard rule).
- Keep the code self-contained (no COMMON_SOURCES dependency), C++23.
- Commit after each verified milestone; never commit with a failing
  Tier-1 gate.

## 15. Deletion (M4 — skeleton only; derive before implementing)

Symmetric setup: `G' = G − e`; cores only fall
(`C_ℓ(G') ⊆ C_ℓ(G)`); `D̄_ℓ = C_ℓ(G) \ C_ℓ(G')` (lost survivors);
faller `x` has `x ∈ D̄_ℓ` exactly for `ℓ ∈ (c'(x), c(x)]`, in particular
`ℓ = c(x)`. The set-version connectivity lemma dualizes cleanly: if
`Z ⊆ D̄_ℓ` and every OLD witness (in `C_ℓ(G)`) of every `z ∈ Z` avoids
`{u,v}` and meets `D̄_ℓ` only inside `Z`, then `C_ℓ(G') ∪ Z` is
ℓ-satisfying in `G'`, contradiction — so fall-chains at level ℓ link
every faller to a destroyed clique (one containing both `u,v`), through
fallers at the SAME level. Consequences to build on: all fallers have
`c <= K_del := min(c(u), c(v))` (the destroyed witness contains `u,v`
inside `C_ℓ(G)` so `ℓ <= K_del`); chain members at level ℓ have
`c >= ℓ`. The discovery test must use PESSIMISM (count only cliques
whose members surely survive) — the exact dual admission/eviction pair
must be derived with the same rigor as §6 (write the lemmas first,
then code), and verified with a deletion harness mirroring §12 (note
deletion reference = remove edge, compare against `core(G−e)` which the
Tier-1 harness already computes). Do NOT improvise shortcuts from the
k-core deletion literature (its ±1 structure does not transfer).

## 16. Scar tissue — the DO-NOT list (each item cost us a real bug)

(a) **No stale-high keys, ever.** A "min" popped with an upper-bound key
    corrupts cores silently (LazyPop UB-key bug; caught only by
    per-vertex diff on ca-HepTh). All keys must be exact at pop time.
(b) **No deferred recounts across distinct levels.** Batch-removing
    pinned across multiple `c`-values and recounting afterwards lets the
    `m`-ratchet skip a vertex's true pop level (produced core 1056 vs
    true 756 on soc-Epinions1 edge (8522,17120)). Batching WITHIN one
    level is safe (Lemma 0 tie-freedom).
(c) **Per-vertex verification only.** Distribution/histogram comparison
    once masked pairwise swaps.
(d) **Strict grow-then-evict** (§6.3). Interleaving breaks Lemma 6.
(e) **Insertion cores never drop; keys never go negative.** Keep both
    guards; treat any firing as a stop-the-line bug.
(f) **Same-level-only optimism is insufficient** (mixed-level rings are
    legal); the static TS-test exists precisely to handle this. Do not
    "simplify" it back to same-level.

## 17. Milestones

- **M1**: bitset counting engine, bit-identical on the 1188 harness
  against the v1 binary (same CHANGED lines, same STATS counts modulo
  timing fields).
- **M2**: Phase 1 + Phase 2 + single Phase-3 peel replacing v1's
  round loop. Gate: Tier-1 1188/1188; report `|C|` vs v1 region sizes
  and `|C|/|R*|`.
- **M3**: measure vs §13 targets; only if p90 misses: leveled-cap
  tightening (§6.5) and/or engine profiling. Re-gate Tier-1 after ANY
  change.
- **M4**: deletion (derive §15 fully first — lemmas before code).
- **M5**: batch insertion (union of seeds, one closure over all new
  edges' chains, one peel — the lemmas generalize: chains end at ANY
  inserted edge's clique) + Tier-2 pokec + larger graphs.

---

# v3 ADDENDUM — killing the discovery flood (post-mortem driven)

**Status**: authoritative extension after v2 measurement. v2's Phases 1/3 are
validated (0/1188 mismatches; Phase-3 peel = 646 µs on the diagnosed edge).
v2's Phase 2 floods: |C|=3064 admitted for a true region of 16 on
soc-Epinions1 s=5 edge (1283,2927) (p2 = 3.4 s, of which eviction 2.99 s);
on s=3 configs the closure exceeds the cap on 78–82% of edges (fallbacks).
Root cause, measured: the static TS-disjunct is vacuous at small thresholds
(ℓ_y = 2..4 admits the whole low-core sea) and on level-homogeneous shells
(everyone's TS exceeds ℓ). The fix below removes the flood AT THE SOURCE
with a new, provably complete level filter, plus one engineering upgrade.

## 18. Active-level theory (the new lever)

Recall `D_ℓ = C_ℓ(G')\C_ℓ(G)`, `Λ = {ℓ : D_ℓ ≠ ∅}` (the ACTIVE levels),
seeds `S = {u,v} ∪ W` (optionally restricted to W-members participating in
at least one new s-clique; for s = 3 that is all of W).

**Lemma 7 (boundary monotonicity).** For a fixed region R, if `b <= b'`
pointwise on `V\R`, then `ccore_{R,b}(x) <= ccore_{R,b'}(x)` for all
`x ∈ R`.
*Proof.* For every ℓ the boundary sets nest:
`{w ∉ R : b(w) >= ℓ} ⊆ {w ∉ R : b'(w) >= ℓ}`, so any valid ℓ-set `D` for
`(R,b)` is valid for `(R,b')`; maximum sets therefore nest. ∎

**Lemma 8 (∞-boundary upper bound).** For ANY region R and `x ∈ R`:
`c'(x) <= ccore_{R,∞}(x)`, where `∞` denotes `b ≡ +∞`.
*Proof.* Let `ℓ = c'(x)` and `D = C_ℓ(G') ∩ R`. Every `z ∈ D` has `>= ℓ`
witnesses inside `C_ℓ(G') ⊆ D ∪ (V\R)`, and with `b ≡ ∞` the level-ℓ
boundary set is all of `V\R`. So D is a valid ℓ-set for `(R,∞)`, and
`x ∈ D`. ∎

**Lemma 9 (seed cut / active-level filter).**
`Λ ⊆ ∪_{x ∈ S} ( c(x), c'(x) ]`, and consequently, with any per-seed upper
bounds `UB(x) >= c'(x)`,
`Λ ⊆ Λ̂ := ∪_{x ∈ S} ( c(x), UB(x) ]`.
Moreover every riser `y ∈ R*` satisfies `c(y)+1 ∈ Λ ⊆ Λ̂`, and every
co-riser member `z ∈ R*` of any witness likewise satisfies
`c(z)+1 ∈ Λ̂`.
*Proof.* `ℓ ∈ Λ ⟹ D_ℓ ≠ ∅`; by Corollary 3a some chain-terminal
`x_m ∈ D_ℓ` has a witness `K ∋ u,v` with `K ⊆ C_ℓ(G')`; `x_m ∈ K` forces
`x_m ∈ {u,v} ∪ W = S`. `x_m ∈ D_ℓ` means `c(x_m) < ℓ <= c'(x_m)`, i.e.
`ℓ ∈ (c(x_m), c'(x_m)]`. The rest: `y ∈ R* ⟹ y ∈ D_{c(y)+1}` ⟹
`c(y)+1 ∈ Λ`; and `c' <= UB` widens intervals only. ∎

**Corollary 9a (second early exit).** If `Λ̂ = ∅` then `R* = ∅`: no core
changes at all; emit STATS and exit (before closure and peel).

**Computing UB (the seed mini-peel).** Run the pinned peel of §4 on region
S with boundary values `b ≡ +∞` — i.e. the boundary never exits; keys drop
only via region pops. By Lemma 1 this computes `ccore_{S,∞}`, which
by Lemma 8 upper-bounds `c'` on S. Capping (hub seeds' exact supports are
the expensive objects we must avoid): compute each seed key with the capped
counter at `K̂_x = c(x) + 256`. A seed whose key saturates is treated as
UB = ∞ (open-topped interval) AND remains permanently alive in the
mini-peel (equivalently it is moved to the ∞-boundary; by Lemma 7 this can
only raise the other seeds' computed values — still sound). Unsaturated
keys may be recomputed from scratch (capped) after each pop; |S| is tiny.
Deltas must never be applied to a saturated key.

## 19. The Λ̂ filter — three hooks into Phase 2 (each completeness-preserving)

Compute Λ̂ (sorted disjoint intervals) BEFORE the closure. Then:

1. **Trigger filter**: skip testing y (and mark it tested — the test is
   static) unless `c(y)+1 ∈ Λ̂`. [Complete by Lemma 9: risers pass.]
2. **PASS member filter**: member z of a counted witness must satisfy
   `c(z) >= ℓ_y` OR (`c(z)+1 ∈ Λ̂` AND `TS(z) >= ℓ_y`). [A member needed
   under optimism is a co-riser, hence in-band by Lemma 9.]
3. **EOS member filter**: z counts via the optimistic branch only if
   `z ∈ C` AND `c(z)+1 ∈ Λ̂` AND `TS(z) >= thr`. [Same argument.]

Seeds u,v stay unconditional members of C. Everything downstream
(eviction order, Phase 3) is unchanged; Lemma 2's hypothesis is still
delivered by (filtered closure ⊇ R*, Theorem 5 + Lemma 9) + (eviction
safety, Lemma 6).

Expected effect (verify empirically, §21): on edge (1283,2927)
(c(u)=766, c(v)=48) every candidate with c below the seeds' interval
bottoms is excluded up front — the measured 3064-flood's low-core mass
dies instantly; on dblp s=3, edges whose seeds do not rise (or rise in
narrow bands) never start a flood, and Λ̂ = ∅ edges exit in Phase 1.5.

## 20. Delta-maintained eviction (Stage B — only if still needed)

If in-band floods remain after §19, replace the recompute-per-recheck
eviction with a delta-maintained cascade using the SAME machinery as the
peel: store ev[x] exact-below-cap (cap = c(x)+1+margin); on each eviction
of z, for each C-neighbor x with ev[x] exact: ev[x] -= (pair count through
{x,z} under x's EOS member-predicate, computed BEFORE removing z's C
membership); saturated ev[x] is recomputed (capped) on first touch.
Decisions (`ev[x] < c(x)+1`) must always be made on exact-at-decision-time
values (§16a analog). Worklist otherwise unchanged.

## 21. Staged execution & measurement plan (measure before engine)

- **Stage A**: implement §18–§19 only. Measure: (i) the four diagnosed
  slow edges of soc-Epinions1 s=5 — (1283,2927), (1749,1822),
  (40200,40202), (148,353): report |C| before/after, p2_us, insert_us;
  (ii) 20-edge smokes on Epinions s=3/s=5; (iii) 300-edge dblp s=3 —
  report fallback count (was 246/300). Decision gate: if p90 collapses
  and fallbacks drop to ~single digits, skip Stage B.
- **Stage B**: §20 if needed. Re-measure.
- **Stage C**: full Tier-1 sweep (all four configs), then the §13
  peel-only acceptance table. Correctness gate remains 0 mismatches.

---

# v4 ADDENDUM — dynamically maintained succinct clique index (index-backed mode)

**Motivation** (user directive + measurement): v3's residual costs are
dominated by SCT recursion inside discovery tests (dblp s=3: ~30k tests,
median 51.7 ms vs 32.6 ms peel-only). The index-backed mode replaces
recursion-based counting with leaf-list evaluation over a MAINTAINED
CPI. Novelty check (2026-07-03): dynamic maximal-clique-set maintenance
exists; batch-dynamic clique COUNT maintenance exists (numbers only);
nobody maintains the succinct index itself. Open lane.

## 22. Additivity of the index under insertion

**Lemma 10 (CPI additivity).** Let F be any exact clique forest for G
(every clique of G encoded at exactly one leaf as H ∪ σ, σ ⊆ Π). Let
T_e be an SCT built on G[W] with held pair {u,v} prepended to every
leaf's H. Then F ⊎ T_e is an exact clique forest for G' = G + e.
*Proof.* Cliques of G' split disjointly into (a) cliques of G (none
contains both u and v since e ∉ G) — encoded exactly once by F, never
by T_e (every T_e clique contains u and v); (b) cliques containing both
u and v — each equals {u,v} ∪ K' with K' a clique of G[W] (members
must be common neighbors), encoded exactly once by T_e (SCT property on
G[W]), never by F. ∎

Insert procedure: build T_e on G[W] (the Phase-1 local SCT, now kept),
append its leaves, extend the V→L incidence lists of {u,v} ∪ (T_e
members), and update the persistent per-vertex total supports via T_e's
per-leaf attribution (held: C(|Π|, k−|H|); pivot: C(|Π|−1, k−|H|−1)).

## 23. Deletion surgery (leaf-local)

Deleting e = (u,v) kills exactly the cliques containing both u,v. For
each leaf L = (H, Π) with both u,v ∈ H ∪ Π (locate via
leaves(u) ∩ leaves(v) using the V→L index), apply exactly one case:
1. u,v ∈ H: delete L (all its cliques contain H ⊇ {u,v}).
2. u ∈ H, v ∈ Π (or symmetric): replace by (H, Π∖{v})
   (survivors are exactly σ ⊆ Π∖{v}).
3. u,v ∈ Π: replace by TWO leaves (H, Π∖{u}) and (H∪{u}, Π∖{u,v})
   (subsets of Π lacking u, plus subsets containing u but not v —
   disjoint and covering all survivors exactly once).
Leaves with at most one of u,v are untouched. Each edit is O(|leaf|)
including incidence-list fixes; total cost O(Σ_{L ∈ leaves(u)∩leaves(v)}
|L|). *Proof of exactness*: case analysis above partitions each leaf's
clique set into dead (⊇{u,v}) and alive, re-encoding the alive part
exactly once. ∎

## 24. Compression drift and amortized rebuild

After T updates the forest is base + edge-trees + surgery fragments;
size Σ_cur drifts above a fresh build's Σ_fresh. Policy (scapegoat
style): track Σ_cur and an estimate of Σ_fresh (e.g. from the last full
build, scaled); when Σ_cur > (1+ε)·Σ_fresh (ε = 0.5 default), rebuild
from scratch (cost = one static build, amortized over the ≥ ε-fraction
of growth updates). Per-root subtree rebuilds are a refinement; measure
first.

## 25. Index-backed counting interface (what discovery/peel consume)

- `TS(z)`: O(1) — maintained persistent support[] (updated in §22/§23).
- Filtered counts (OS/EOS/keys; member predicate `pred`):
  Σ over L ∈ leaves(x) with H(L) all-pass of C(#{π ∈ Π(L) pass}, s−|H|).
  Per-leaf pivot scan is O(|Π|); for threshold predicates (c(z) >= thr)
  keep each leaf's Π sorted by core value → binary search per leaf.
  (Core values change only inside R* per update — resort only touched
  leaves.)
- Peel with alive-filtering: V3's per-leaf counter machinery
  (remainPivots/needPivot) scoped to leaves(region) — index-native
  scoped peel.
Persistent state (index-backed mode): core[] + adjacency + forest +
V→L incidences + support[]. Memory O(Σ) — the index-free mode (v3)
remains the light alternative; the paper presents BOTH backends over
the same discovery theory (Λ̂/closure/eviction/pinned peel unchanged).

## 26. Kill-or-confirm experiments (build first, decide from data)

- **E1 (additivity, bit-exact)**: standalone driver: per sampled edge e
  and s ∈ {3,5}: per-vertex s-clique counts of G' computed (a) fresh
  SCT on G' vs (b) SCT on G plus T_e attribution. Require exact
  equality on every vertex, ≥100 edges × 2 graphs. Kills Lemma 10's
  implementation risk.
- **E2 (surgery cost)**: one SDCT walk per graph recording, for the 300
  sampled edges, |leaves(u) ∩ leaves(v)| and Σ|L| over those leaves
  (single pass, callback-based; no per-pair walks). Distribution decides
  deletion viability.
- **E3 (index value pricing)**: from Stage-C phase breakdowns, the SCT
  recursion share of discovery time per config = the exact upper bound
  on what §25 can erase. Compute projected medians/p90s.
Decision rule: E1 exact AND E2 median surgery cost ≲ discovery cost AND
E3 projects peel-only wins on all four configs ⟹ implement v4;
otherwise index-free v3 stays the primary and v4 is scoped down.

## 27. v4 implementation plan (decision: GO, per §26 — E1 exact, E2 median
cheap, E3 96.8-99.8% of residual cost is discovery counting; on
epinions s=5 the peel co-dominates at 48%, so the index must also serve
peel keys/deltas, not just discovery)

Milestones, each gated:
- **V4.1 index at load**: build the level-s clique forest (leaves =
  (H, Π) pairs + V→L incidence CSR) at startup (untimed, like graph
  load), plus maintained support[] (per-vertex totals from leaf
  attribution). Storage: flat arrays, leaf Π sorted by core value.
- **V4.2 index-backed discovery**: TS(z) = support[] O(1); OS/EOS =
  per-leaf evaluation over leaves(y): H∖{y} all-pass check + count of
  qualifying Π members (core-threshold part by binary search on the
  sorted Π; the TS-disjunct part by O(1)-per-member predicate scan).
  Gate: Tier-1 all four configs 0 failures AND CHANGED lines
  bit-identical to the Stage-B binary.
- **V4.3 index-native peel**: keys/deltas served from leaves(region)
  with V3-style per-leaf counters (remainPivots/needPivot) instead of
  per-call universe rebuild + recursion. Same gate.
- **V4.4 index maintenance across updates**: append EdgeTree (Lemma 10)
  after each insert + update support[] and incidences. Gate: Tier-2
  streaming (pokec s=3: remove 50 random edges, insert back one at a
  time REUSING the maintained index; final state must equal core(G)
  per-vertex) — this exercises Lemma 10 end-to-end.
- **V4.5**: final dual-mode table (index-backed vs index-free vs both
  baselines); deletion surgery (§23) after that.
