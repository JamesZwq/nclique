# R=1 Memory-Efficient CPI: Contributions for TODS Extension

> This document contains all theoretical results, algorithm descriptions, complexity analyses,
> and experimental evidence for the R=1 optimizations (ST and ST_V2) over the SIGMOD baseline.
> Intended as source material for the TODS journal paper.

---

## 1. Problem Setting Recap

We consider the $(1,s)$-nucleus decomposition problem: given an undirected graph $G=(V,E)$ and an integer $s \ge 2$, compute for every vertex $v \in V$ its $(1,s)$-core number $\kappa_s(v)$, defined as the maximum value $k$ such that $v$ belongs to a subgraph where every vertex participates in at least $k$ distinct $s$-cliques.

The SIGMOD paper introduces the **Clique Path Index (CPI)**, a compact bipartite representation of the $s$-clique search space as a set of paths $\mathcal{P} = \{P_1, P_2, \ldots, P_L\}$. Each path $P = (V_h(P), V_p(P))$ consists of a **hold set** $V_h(P)$ and a **pivot set** $V_p(P)$, encoding $\binom{|V_p(P)|}{s - |V_h(P)|}$ distinct $s$-cliques via Lemma 3.3 of the SIGMOD paper.

Throughout this document:
- $L$ = number of CPI paths
- $\Sigma = \sum_{P \in \mathcal{P}} |P|$ = total vertex-path incidences (sum of path sizes)
- $n = |V|$, $m = |E|$
- $p(P) = |V_p(P)|$, $h(P) = |V_h(P)|$, $\eta(P) = s - h(P)$ (the "needPivot" count)

---

## 2. The SIGMOD Baseline Algorithm (R=1)

The SIGMOD paper's R=1 implementation (`NCliqueVertexCoreDecomposition`) performs peeling as follows:

### 2.1 Data Structures

| Structure | Type | Space |
|---|---|---|
| `tree.adj_list[L]` | `vector<vector<TreeGraphNode>>` — mutable CPI paths | $O(\Sigma)$, ~8 bytes/entry |
| `treeGraphV.adj_list[v]` | `vector<robin_hood::unordered_flat_set<TreeGraphNode>>` — vertex→path hash index | $O(\Sigma)$, ~80 bytes/entry (hash overhead) |
| `countingV[v]` | `double[]` — per-vertex $s$-support | $O(n)$, 8 bytes/vertex |
| `heap` + `heapHandles[v]` | `boost::d_ary_heap<8>` — mutable min-heap | $O(n)$, ~32 bytes/vertex |
| `leafRmInfo[L]` | `StaticVector<LeafRmInfo>` — per-path removal tracking | $O(L)$, up to ~1600 bytes/path (reserves capacity for removedPivots) |

**Total space**: $O(\Sigma + n)$, but with dominant constant factor **~88·Σ + 40·n bytes**.

### 2.2 Per-Iteration Operations

Each peeling iteration:

1. **Extract-min**: `heap.pop()` — $O(\log n)$
2. **Find affected paths**: For each removed vertex $v$, iterate `treeGraphV.getNbr(v)` — hash set iteration, $O(\deg_{\mathcal{P}}(v))$ but cache-unfriendly
3. **Classify each path**: Check if hold vertex removed (Case A: path dies) or pivot removed (Case C: path shrinks)
4. **Compute support delta**: Using nCr formula (SIGMOD Theorem 5)
5. **Update supports**: For each surviving vertex on affected paths, subtract delta from `countingV[v]`
6. **Update heap**: `heap.update(heapHandles[v])` — $O(\log n)$ per affected vertex
7. **Mutate tree**: `tree.removeNode(leafId)` or `tree.removeNbrs(leafId, removedPivots)` + `treeGraphV.removeNbr(v, leafId)` — hash set erase, $O(1)$ amortized but with poor cache behavior

### 2.3 Total Time Complexity

$$T_{\text{SIGMOD}} = O(\Sigma \cdot \log n)$$

The $\log n$ factor comes from heap operations: each vertex receives $O(\deg_{\mathcal{P}}(v))$ support updates across all iterations, and each update triggers `heap.update()` costing $O(\log n)$.

Summing over all vertices: $\sum_v \deg_{\mathcal{P}}(v) = \Sigma$, giving total heap update cost $O(\Sigma \log n)$.

---

## 3. Contribution 1: Immutable-CPI Peeling (ST)

### 3.1 Key Observation

In the SIGMOD algorithm, tree mutation (Step 7) is the most complex operation: it involves dynamic memory reallocation (`removeNode` deallocates a vector, `removeNbrs` shifts elements), hash set erasure (`treeGraphV.removeNbr` probes and deletes from a hash table), and bookkeeping (`leafRmInfo` records which pivots were removed, requiring sort + dedup).

**We observe that for R=1, none of this mutation is necessary.** The support delta for every vertex on an affected path can be computed from just two integers — the pivot count before and after the batch removal — without modifying the CPI structure.

### 3.2 Formal Results

**Lemma 1 (Counter-Based Delta for R=1).**
Let $P = (V_h(P), V_p(P))$ be a CPI path with $p = |V_p(P)|$ pivots and $\eta = s - |V_h(P)|$ needPivot.
Suppose $d$ pivots are removed from $P$ in a batch (and no hold vertex is removed). Then the change in $s$-support for each remaining vertex is:

$$\Delta_{\text{hold}} = \binom{p}{\eta} - \binom{p - d}{\eta}$$

$$\Delta_{\text{pivot}} = \binom{p-1}{\eta-1} - \binom{p-d-1}{\eta-1}$$

where $\Delta_{\text{hold}}$ applies to every hold vertex and $\Delta_{\text{pivot}}$ to every surviving pivot vertex.

*Proof.* By the SIGMOD paper's Lemma 3.3, each hold vertex $v \in V_h(P)$ participates in all $\binom{p}{\eta}$ encoded $s$-cliques before removal. After removing $d$ pivots, the remaining pivot count is $p' = p - d$, and the new count is $\binom{p'}{\eta}$. The delta is their difference: $\Delta_{\text{hold}} = \binom{p}{\eta} - \binom{p'}{\eta}$.

For a surviving pivot vertex $u \in V_p(P) \setminus V_{\text{rm}}$: before removal, $u$ participates in $\binom{p-1}{\eta-1}$ of the $\binom{p}{\eta}$ cliques (those that include $u$ and choose $\eta - 1$ other pivots from the remaining $p - 1$). After removal, the count becomes $\binom{p'-1}{\eta-1} = \binom{p-d-1}{\eta-1}$. □

**Lemma 2 (Counter-Based Delta for Path Death).**
If a hold vertex is removed from $P$, or if $d$ pivots are removed such that $p - d < \eta$ (insufficient pivots), the path dies. The support subtracted from each remaining vertex is the full contribution:

$$\Delta_{\text{hold}}^{\dagger} = \binom{p}{\eta}, \qquad \Delta_{\text{pivot}}^{\dagger} = \binom{p-1}{\eta-1}$$

*Proof.* By Theorem 5 of the SIGMOD paper (Vertex Removal), removing a hold vertex eliminates all $s$-cliques on the path. The full contribution of $P$ to each hold vertex is $\binom{p}{\eta}$, and to each pivot vertex is $\binom{p-1}{\eta-1}$. These are subtracted in full. □

**Theorem 1 (Immutable-CPI Correctness for R=1).**
For $(1,s)$-nucleus decomposition, the CPI structure $(V_h(P), V_p(P))$ for every path $P$ remains unchanged throughout the entire peeling process. The core number $\kappa_s(v)$ computed using only integer counters $(p_P, \eta_P, \text{alive}_P)$ per path equals the core number computed by the SIGMOD algorithm with tree mutation.

*Proof.* We show that tree mutation in the SIGMOD algorithm is redundant for R=1: it changes the tree state but does not affect the computed core numbers.

**Claim**: At any point during peeling, the support of an active vertex $v$ can be computed as:
$$\deg_s(v) = \sum_{\substack{P \in \mathcal{P}: v \in P \\ P \text{ alive}}} w(v, P)$$
where $w(v, P) = \binom{p_P^{\text{current}}}{\eta_P}$ if $v \in V_h(P)$ and $w(v, P) = \binom{p_P^{\text{current}} - 1}{\eta_P - 1}$ if $v \in V_p(P)$, and $p_P^{\text{current}}$ is the current remaining pivot count maintained by the counter.

**Why the tree is not needed**: In R=1, peeling removes **vertices** (1-cliques). By Theorem 5 of the SIGMOD paper:
- Case A (hold vertex removed): Path dies → alive$_P$ ← false. No structural change needed.
- Case C (pivot vertex removed): Pivot removed from path → $p_P$ ← $p_P - 1$. The hold set $V_h(P)$ is unchanged, and the *identity* of remaining pivots does not matter for the nCr formula — only their *count* matters.

**Key insight**: The nCr delta formulas (Lemmas 1 and 2) depend only on $(p_P^{\text{before}}, p_P^{\text{after}}, \eta_P)$, not on *which* pivots were removed. Therefore, maintaining $(p_P, \eta_P, \text{alive}_P)$ as integer counters suffices.

**No Case B exists for R=1**: Vertex removal never creates new paths (unlike edge removal for R=2, which has the pivot-pivot split case). The set of paths only shrinks.

**Formal equivalence**: Let $\kappa_s^{\text{mut}}(v)$ be the core number from the SIGMOD mutable-tree algorithm, and $\kappa_s^{\text{imm}}(v)$ be our immutable-counter algorithm. We prove $\kappa_s^{\text{mut}}(v) = \kappa_s^{\text{imm}}(v)$ by induction on peeling iterations.

*Base case*: Initially, both algorithms compute the same support $\deg_s(v)$ from the same CPI paths.

*Inductive step*: Suppose after $t$ iterations, both algorithms have the same set of active vertices and the same support values. In iteration $t+1$:
- Both pop the same minimum-support vertex $v^*$ (same bucket/heap state).
- SIGMOD mutates tree: removeNode or removeNbrs. Our algorithm: decrement $p_P$ or set alive$_P$ = false.
- By Lemmas 1–2, the support deltas computed by both algorithms are identical.
- Therefore, after the update, both algorithms have the same support values.

By induction, all core numbers are identical. □

**Corollary 1 (No Dynamic Memory Operations).**
The immutable-CPI algorithm performs zero dynamic memory allocations, zero hash table operations, and zero vector resize/shift operations during the entire peeling phase. All data structures are read-only arrays with integer counter updates.

### 3.3 Algorithmic Consequences

**Theorem 2 (Heap Elimination).**
When the CPI is immutable, the mutable min-heap (boost::d_ary_heap with O(log n) decrease-key) can be replaced by a bucket array with O(1) amortized bucket-move operations.

*Proof.* The SIGMOD algorithm requires a mutable heap because `heap.update(handle)` must be called after each support decrease. With immutable CPI, support values are non-negative integers (they are sums of nCr values, which are integers). We use a standard bucket decomposition: `buckets[k]` stores all vertices with current support $k$. After a support decrease from $k$ to $k'$:
- Remove vertex from `buckets[k]` in O(1) (swap with last element + pop).
- Insert into `buckets[k']` in O(1) (append).

This replaces O(log n) heap operations with O(1) bucket moves. □

**Theorem 3 (Hash-Free Vertex-Path Index).**
When the CPI is immutable, the vertex→path index (`treeGraphV`) can be implemented as a static CSR (Compressed Sparse Row) array instead of a per-vertex hash set, reducing per-entry space from ~80 bytes to 5 bytes.

*Proof.* The SIGMOD algorithm uses `robin_hood::unordered_flat_set<TreeGraphNode>` for `treeGraphV` because it needs to support:
1. `removeNbr(v, leafId)` — delete a path from v's neighbor set
2. `addNbr(v, leafId)` — add a path (for Case B, R≥2 only)

For R=1: (1) is unnecessary because the CPI is immutable (no path removal from the index); (2) is impossible (no new paths created). Therefore, `treeGraphV` only needs to support **iteration** over v's paths, which is exactly what a static CSR provides.

**CSR layout**:
- `vtxLeafOff[v]` to `vtxLeafOff[v+1]`: range of entries for vertex $v$ in the data array
- `vtxLeafData[i]` = `{leafId: uint32, isPivot: uint8}` = 5 bytes per entry

**Hash set layout** (robin_hood): each entry occupies ~80 bytes including:
- 8 bytes for the `TreeGraphNode` key
- ~72 bytes overhead (load factor ~87%, metadata byte per bucket, alignment padding, power-of-2 bucket count with empty slots)

Space reduction: **80 / 5 = 16× per entry**. □

### 3.4 Complexity Summary (ST)

| | SIGMOD Baseline | ST (Immutable-CPI) |
|---|---|---|
| **Init (support computation)** | $O(\Sigma)$ | $O(\Sigma)$ |
| **Init (vertex→path index)** | $O(\Sigma)$ hash inserts | $O(\Sigma)$ CSR build (2 linear scans) |
| **Init (priority queue)** | $O(n \log n)$ heap build | $O(n)$ bucket init |
| **Peeling (per vertex update)** | $O(\log n)$ heap decrease-key | $O(1)$ bucket move |
| **Peeling (tree mutation)** | $O(1)$ amortized hash erase | **Eliminated** |
| **Peeling total** | $O(\Sigma \cdot \log n)$ | $O(\Sigma)$ |
| **Total time** | $O(\Sigma \cdot \log n)$ | $O(\Sigma)$ |
| **Space (vertex→path)** | ~80·Σ bytes | ~5·Σ bytes |
| **Space (per-path metadata)** | ~1600·L bytes (leafRmInfo) | 9·L bytes (3 integers) |
| **Space (priority queue)** | ~32·n bytes (heap + handles) | ~12·n bytes (bucket arrays) |

---

## 4. Contribution 2: Tree-Free CPI Construction (ST_V2)

### 4.1 Motivation

ST (Section 3) makes the CPI immutable during peeling but still materializes the tree during construction. The tree `DynamicGraph<TreeGraphNode>.adj_list` stores every path as a `vector<TreeGraphNode>`, consuming 8 bytes per vertex-path entry (the `TreeGraphNode` struct: 63-bit vertex ID + 1-bit isPivot flag, padded to 8 bytes).

After the tree is built, ST constructs the flat CSR index from it (one additional O(Σ) pass), then computes initial supports (another O(Σ) pass). The tree itself is never accessed again during peeling — it's only used as an intermediate representation.

**Observation**: If we can build the CSR index and compute supports *during* CPI construction (inside the SDCT Bron-Kerbosch recursion), the tree never needs to exist.

### 4.2 Construction-Integrated Indexing

**Definition (Augmented SDCT).** An Augmented SDCT is a modified Succinct Degeneracy Clique Tree enumeration that, instead of storing each leaf as a tree node, invokes a callback function $\texttt{onLeaf}(\text{leafId}, V_h, V_p)$ at each BK-recursion leaf. The callback receives the leaf's hold and pivot vertex sets and a globally unique leaf identifier.

**Theorem 4 (Tree-Free CPI Construction).**
The CPI tree, the vertex→path CSR index, the path→vertex CSR index, and the initial per-vertex $s$-support can all be constructed in a single SDCT pass without materializing the tree, using $O(\Sigma)$ total time and $O(\Sigma + n + L)$ space.

*Proof.* We modify the SDCT algorithm to invoke a callback at each leaf instead of storing the leaf. The callback:

1. Records `(vertex, leafId, isPivot)` triples into a COO (Coordinate) buffer — O(1) per vertex-path entry.
2. Computes initial support for each vertex inline:
   - Hold vertex $v$: $\deg_s(v) \mathrel{+}= \binom{|V_p|}{\eta}$
   - Pivot vertex $v$: $\deg_s(v) \mathrel{+}= \binom{|V_p|-1}{\eta-1}$
3. Records per-leaf metadata: $p_P = |V_p|$, $\eta_P = s - |V_h|$.

After the SDCT completes, the COO buffer is converted to dual CSR in O(Σ) time (count + prefix-sum + scatter, two passes).

**Three passes eliminated**:
- Pass 1 (tree → treeGraphV construction): hash set insertions — **eliminated** (CSR built from COO)
- Pass 2 (tree → countingPerVertex): support scan — **eliminated** (inline in callback)
- Pass 3 (tree → flat CSR): CSR build from tree — **eliminated** (CSR built from COO)

The tree itself is never stored.

**Implementation detail (AugCtx)**: The BK recursion can reach depths of hundreds of levels. Using `std::function` for the callback would add ~40 bytes of stack overhead per recursion level, risking stack overflow on large graphs. We use a function-pointer trampoline (`AugCtx` context struct) that provides zero-overhead type erasure: the callback is stored as a raw function pointer + void* context, adding only 16 bytes per frame. □

### 4.3 Space Analysis

| Component | SIGMOD | ST | ST_V2 |
|---|---|---|---|
| CPI tree (`adj_list`) | 8·Σ bytes | 8·Σ bytes | **0 bytes** |
| Vertex→path index | ~80·Σ bytes (hash set) | 5·Σ + 4n bytes (CSR) | 5·Σ + 4n bytes (CSR) |
| Path→vertex index | — (implicitly via tree) | — (implicitly via tree) | 5·Σ + 4L bytes (CSR) |
| Per-path metadata | ~1600·L bytes | 9·L bytes | 9·L bytes |
| Per-vertex support | 8n bytes | 8n bytes | 8n bytes |
| Priority queue | ~32n bytes | ~12n bytes | ~12n bytes |
| **Total dominant term** | **~88·Σ** | **~13·Σ** | **~10·Σ** |
| **Reduction factor** | 1× | **6.8×** | **8.8×** |

**Remark**: The space reduction factor on the Σ-term ranges from 6.8× (ST, tree retained but hash eliminated) to 8.8× (ST_V2, tree-free). In practice, Σ dominates for large $s$: for the `facebook_combined` graph with $s=5$, Σ ≈ 14.4M entries, and the space difference is 88×14.4M ≈ 1.27 GB (SIGMOD) vs 10×14.4M ≈ 144 MB (ST_V2).

### 4.4 Time Analysis

| Phase | SIGMOD | ST | ST_V2 |
|---|---|---|---|
| SDCT enumeration | $O(\Sigma)$ | $O(\Sigma)$ | $O(\Sigma)$ |
| Tree materialization | $O(\Sigma)$ | $O(\Sigma)$ | **0** |
| treeGraphV construction | $O(\Sigma)$ hash inserts | $O(\Sigma)$ CSR build | **0** (fused) |
| countingPerVertex | $O(\Sigma)$ scan | $O(\Sigma)$ scan | **0** (fused) |
| CSR build from tree | — | $O(\Sigma)$ | **0** (from COO) |
| COO → dual CSR | — | — | $O(\Sigma)$ |
| **Init passes (post-SDCT)** | **3 passes** | **2 passes** | **1 pass** |
| Peeling | $O(\Sigma \log n)$ | $O(\Sigma)$ | $O(\Sigma)$ |
| **Total** | $O(\Sigma \log n)$ | $O(\Sigma)$ | $O(\Sigma)$ |

**Remark on SDCT time**: The SDCT enumeration itself takes $O(\Sigma)$ in all three algorithms and dominates the initialization cost. For large graphs, SDCT accounts for 70–90% of total runtime (see Section 6). The peeling phase is where ST/ST_V2 achieve asymptotic improvement over SIGMOD.

---

## 5. Per-Phase Peeling Algorithm

We present the complete peeling algorithm shared by ST and ST_V2. The only difference is the data source: ST reads from `tree.adj_list[leafId]` in Phase 2, while ST_V2 reads from `leafVtxData[leafVtxOff[leafId]..leafVtxOff[leafId+1]]`.

### Algorithm: Immutable-CPI Peeling for (1,s)-Nucleus Decomposition

**Input**: Dual CSR (vtxLeafOff/Data, leafVtxOff/Data), per-path metadata (leafRemainPivots, leafNeedPivot), initial support (countingV).

**Output**: Core number array coreV[v] for each vertex.

```
Initialize:
  For each path P: leafRemainPivots[P] ← |V_p(P)|
                    leafNeedPivot[P] ← s - |V_h(P)|
                    leafAlive[P] ← (leafNeedPivot[P] ≤ leafRemainPivots[P])
  For each vertex v with countingV[v] > 0: insert v into buckets[countingV[v]]
  curLevel ← 0

While active vertices remain:
  // --- Batch Extract ---
  Advance curLevel to first non-empty bucket
  Pop all vertices v with countingV[v] ≤ curLevel → batch V_rm
  Set coreV[v] ← curLevel for each v ∈ V_rm

  // --- Phase 1: Find Affected Paths (via Vertex→Path CSR) ---
  For each v ∈ V_rm:
    For each (pathId, isPivot) in vtxLeafData[vtxLeafOff[v]..vtxLeafOff[v+1]]:
      If not leafAlive[pathId]: skip
      Mark pathId as affected
      If isPivot = false: mark path as dying (Case A)
      Else: leafRemovedPivots[pathId]++

  // --- Phase 2: Compute Deltas (via Path→Vertex CSR) ---
  For each affected path P:
    old_rp ← leafRemainPivots[P]
    η ← leafNeedPivot[P]
    d ← leafRemovedPivots[P]
    new_rp ← old_rp - d
    dies ← (path marked dying) OR (new_rp < η)

    If dies:  // Lemma 2
      Δ_hold ← C(old_rp, η)
      Δ_pivot ← C(old_rp - 1, η - 1)
    Else:  // Lemma 1
      Δ_hold ← C(old_rp, η) - C(new_rp, η)
      Δ_pivot ← C(old_rp - 1, η - 1) - C(new_rp - 1, η - 1)

    For each (vertex, isPivot) in leafVtxData[leafVtxOff[P]..leafVtxOff[P+1]]:
      If vertex not active: skip
      δ ← isPivot ? Δ_pivot : Δ_hold
      countingV[vertex] -= δ
      Mark vertex as dirty

    leafRemainPivots[P] ← new_rp
    If dies: leafAlive[P] ← false

  // --- Phase 3: Batch Bucket Moves ---
  For each dirty vertex v:
    Move v from buckets[old_bucket] to buckets[countingV[v]]  // O(1) each

  Clear batch tracking arrays
```

### Correctness

Follows from Theorem 1 (Immutable-CPI Correctness). Each iteration produces the same set of extracted vertices and the same support updates as the SIGMOD mutable-tree algorithm.

### Time per iteration

Let $B$ = batch size, $A$ = number of affected paths, $\Sigma_A$ = total vertex-path entries on affected paths.

- Phase 1: $O(\sum_{v \in B} \deg_{\mathcal{P}}(v))$ — CSR scan, sequential access
- Phase 2: $O(\Sigma_A)$ — CSR scan, sequential access
- Phase 3: $O(|\text{dirty vertices}|)$ — O(1) per bucket move

Total across all iterations: $O(\Sigma)$ (each vertex-path entry is touched at most twice — once in Phase 1 when the vertex is removed, once in Phase 2 when the path is affected).

---

## 6. Experimental Results

### 6.1 Setup

- **Machine**: Single-threaded, OMP_NUM_THREADS=1
- **Build**: C++23, Clang/G++, -O3
- **Production mode**: No `PIVOTER_COMPARE`, no reference tree for ST_V2
- **Measurement**: Wall-clock time via `/usr/bin/time -l`, peak RSS from OS

### 6.2 End-to-End Wall Time (ms)

| Graph | n | s | SIGMOD | ST | ST_V2 | Speedup (V2/SIGMOD) |
|---|---|---|---|---|---|---|
| email-Eu-core | 1,005 | 3 | 60 | — | 40 | 1.5× |
| | | 4 | 70 | — | 40 | 1.8× |
| | | 5 | 70 | — | 50 | 1.4× |
| facebook | 4,039 | 3 | 5,670 | — | 2,930 | **1.9×** |
| | | 4 | 20,250 | — | 10,410 | **1.9×** |
| | | 5 | 63,950 | — | 34,720 | **1.8×** |
| com-dblp | 317,081 | 3 | 1,030 | — | 790 | 1.3× |
| | | 4 | 1,020 | — | 740 | 1.4× |
| | | 5 | 1,030 | — | 830 | 1.2× |
| com-youtube | 1,134,890 | 3 | 6,820 | — | 3,770 | 1.8× |
| | | 4 | 6,960 | — | 4,030 | 1.7× |
| | | 5 | 6,720 | — | 4,120 | 1.6× |
| web-Stanford | 281,903 | 3 | 9,400 | — | 3,870 | **2.4×** |
| | | 4 | 8,720 | — | 3,620 | **2.4×** |
| | | 5 | 8,490 | — | 3,810 | **2.2×** |

### 6.3 Peak RSS Memory (MB)

| Graph | s | SIGMOD | ST_V2 | Reduction |
|---|---|---|---|---|
| facebook | 3 | 133 | 73 | **1.8×** |
| | 4 | 962 | 623 | **1.5×** |
| | 5 | 5,670 | 4,097 | **1.4×** |
| com-dblp | 3 | 112 | 101 | 1.1× |
| | 4 | 112 | 80 | 1.4× |
| | 5 | 258 | 226 | 1.1× |
| com-youtube | 3 | 533 | 373 | 1.4× |
| | 4 | 628 | 399 | **1.6×** |
| | 5 | 665 | 372 | **1.8×** |
| web-Stanford | 3 | 364 | 303 | 1.2× |
| | 4 | 464 | 415 | 1.1× |
| | 5 | 1,042 | 989 | 1.1× |

### 6.4 Peeling-Only Speedup

The peeling phase (excluding SDCT construction) shows the clearest speedup, as it isolates the effect of immutable-CPI + bucket array vs mutable-tree + heap:

| Graph | s=3 | s=4 | s=5 |
|---|---|---|---|
| facebook | 3.2× | 4.7× | **5.2×** |
| com-dblp | 2.1× | 3.1× | 2.5× |
| com-youtube | 2.7× | 2.8× | 2.5× |
| web-Stanford | 2.8× | 2.3× | 1.7× |

### 6.5 Breakdown Analysis

The SDCT construction time is nearly identical between SIGMOD and ST_V2 (~2–3 seconds on large graphs). The wall-time speedup comes from two sources:

1. **SIGMOD builds `treeGraphV` (hash set) + runs SDCT_Par7 verification** — ST_V2 skips both since `needsSDCT()` returns false. On web-Stanford this alone saves ~5 seconds.
2. **Peeling is 2–5× faster**: O(1) bucket moves vs O(log n) heap updates, sequential CSR access vs random hash probes.

---

## 7. Summary of Theoretical Results

| # | Result | Type | Key Implication |
|---|---|---|---|
| Lemma 1 | Counter-based support delta (pivot removal) | Analytical | nCr difference computable from 2 integers |
| Lemma 2 | Counter-based support delta (path death) | Analytical | Full contribution computable from 2 integers |
| **Theorem 1** | **Immutable-CPI correctness for R=1** | **Algorithmic** | **CPI needs no mutation during peeling** |
| Corollary 1 | Zero dynamic memory during peeling | Algorithmic | No alloc/free/hash ops |
| Theorem 2 | Heap → bucket array | Complexity | O(Σ log n) → O(Σ) peeling time |
| Theorem 3 | Hash set → static CSR | Space | ~80 → 5 bytes per entry |
| **Theorem 4** | **Tree-free CPI construction** | **Systemic** | **Tree never materialized, 3 passes → 1 pass** |

### Why this does not generalize to R≥2

For R=2, the SIGMOD paper's Theorem 4 (Edge Removal) identifies three cases:
- Case 1 (hold-hold): path dies — analogous to R=1 Case A
- Case 2 (hold-pivot): pivot removed — analogous to R=1 Case C
- **Case 3 (pivot-pivot): path splits into two new paths** — **no R=1 analogue**

Case 3 creates new CPI paths with new hold/pivot structures. Subsequent Case 1/2 operations must reference these new paths. Therefore, the CPI cannot remain immutable: it must support path insertion (for split children) and path deletion (for the split parent). This fundamentally requires mutable data structures.

For R≥3, the general path-splitting theorem (SIGMOD Theorem 3) involves BK-based enumeration of sub-paths, producing an arbitrary number of new paths. The same argument applies.

**Corollary**: Immutable-CPI peeling is an R=1-specific optimization. It exploits the structural property that vertex removal (R=1) only kills or shrinks paths, never creates new ones.
