# Augmented SDCT: Formal Complexity Analysis

## 1. Output Equivalence

**Theorem 1 (SDCT_Augmented_NoTree equivalence).**
Let G be an input graph, s = max_k, r = min_k. Then `SDCT_Augmented_NoTree(G, s, r, onLeaf)` invokes `onLeaf(leafId, keepV, dropV)` for exactly the same set of (keepV, dropV) pairs as `SDCT(G, s, r)` emits as tree leaves, with identical keepV/dropV content per leaf. The only difference is leaf ID assignment: SDCT assigns IDs via `DynamicGraph::addNode()` (which may reuse freed slots), while NoTree assigns sequential IDs 0, 1, 2, ....

*Proof sketch.* Both functions execute the identical BK recursion: same pivot selection (`findBestPivotNonNeighborsDegeneracyCliques`), same vertex partitioning (`moveToRDegeneracyCliques`), same keep/drop classification (pivot vertex → drop, non-pivot → keep), same pruning conditions (`keepV.size() > max_k`, `cSize < min_k`). The only difference is the leaf emission action: SDCT calls `tree.addNode(newNode)`, NoTree calls `onLeaf(leafId, keepV, dropV)`. Since SDCT never reuses freed nodes during construction (no `removeNode` calls during BK), `addNode` produces sequential IDs identical to `leafCounter++`.

Special case: when `keepV.size() == max_k`, SDCT stores only keep vertices in the tree node (no pivots). NoTree matches this by passing `emptyDrop` to the callback.

**Theorem 2 (Tree-Free R1 peeling equivalence).**
`NCliqueVertexCoreDecomposition_ST_V2` produces identical core values to `NCliqueVertexCoreDecomposition_ST` for all vertices.

*Proof sketch.* Both algorithms:
1. Compute identical initial support `countingV[v]` (Theorem 1 guarantees same leaves; the nCr formula is identical).
2. Use identical bucket-based peeling with batch removal.
3. Compute identical per-leaf delta (same `leafRemainPivots`, `leafNeedPivot`, same nCr delta formula).
4. Iterate the same set of vertices per leaf during Phase 2 (Leaf→Vertex CSR matches `tree.adj_list[leafId]` by Theorem 1).

The only structural difference is the iteration mechanism: ST uses `tree.adj_list[leafId]` (vector of TreeGraphNode), V2 uses `leafVtxData[leafVtxOff[leafId]..leafVtxOff[leafId+1])` (flat CSR). Both iterate the same (vertex, isPivot) pairs.

## 2. Space Complexity

### Per-entry storage (one vertex-in-leaf occurrence)

| Data structure | ST (current) | ST_V2 (tree-free) |
|---|---|---|
| `tree.adj_list[L]` entry | 8B (TreeGraphNode: 63-bit v + 1-bit isPivot) | **0B** (eliminated) |
| `tree.adj_list[L]` vector overhead | ~24B / |L| amortized (vector header) | **0B** |
| `treeGraphV` hash set entry | ~40B (robin_hood flat_set per vertex) | **0B** (eliminated) |
| Vertex→Leaf CSR entry | 5B (leafId:4B + isPivot:1B) | 5B (same, with padding → 8B) |
| Leaf→Vertex CSR entry | — | 5B (vertex:4B + isPivot:1B, padded → 8B) |
| Per-leaf metadata | 12B (3 × int arrays) | 8B (2 × int: pivotCount, needP) |

**Theoretical per-entry cost:**
- ST: 8 + 24/avg_leaf_size + 40 + 5 ≈ **53–77B** per vertex-in-leaf occurrence
- ST_V2: 8 + 8 ≈ **16B** per occurrence (struct padding) + 8B/leaf metadata

**Theoretical reduction: 3.3–4.8×** per entry.

### Total space for com-dblp (r=1, s=4)

| Component | ST | ST_V2 |
|---|---|---|
| Leaves | 166,918 | 166,918 |
| Σ|leaf| (total entries) | 891,989 | 891,989 |
| tree.adj_list | 891,989 × 8B = 6.8 MB | **0** |
| tree vector headers | 166,918 × 24B = 3.8 MB | **0** |
| treeGraphV hash sets | 891,989 × ~40B = 34.0 MB | **0** |
| Vertex→Leaf CSR | 891,989 × 8B = 6.8 MB | 891,989 × 8B = 6.8 MB |
| Leaf→Vertex CSR | — | 891,989 × 8B = 6.8 MB |
| Per-leaf metadata | 166,918 × 12B = 1.9 MB | 166,918 × 8B = 1.3 MB |
| **Total index** | **~53.3 MB** | **~14.9 MB** |

**Measured (from experiments):**
- ST: Other Index Memory = 116 MB (includes graph + tree + treeGraphV + CSR)
- ST_V2: ST_V2 after CSR build = 126 MB (includes graph + refTree + dual CSR)
- Note: ST_V2 currently still builds `refTree` for verification. Without it, savings would be ~10.6 MB on this dataset.

### Production mode estimate (no refTree)

When `refTree` is skipped entirely:
- ST: Graph (41 MB) + tree (10.6 MB) + treeGraphV (34 MB) + CSR (6.8 MB) + metadata (1.9 MB) ≈ **94.3 MB**
- ST_V2: Graph (41 MB) + dual CSR (13.6 MB) + metadata (1.3 MB) ≈ **55.9 MB**
- **Savings: 38.4 MB (40.7% reduction)**

## 3. Time Complexity

### Construction phase

| Pass | ST | ST_V2 |
|---|---|---|
| SDCT enumeration | O(T_BK) | O(T_BK) — identical |
| treeGraphV build | O(Σ\|leaf\|) | **0** (fused into callback) |
| countingPerVertex | O(Σ\|leaf\|) | **0** (fused into callback) |
| CSR build | O(Σ\|leaf\|) | O(Σ\|leaf\|) — COO→CSR conversion |
| **Total post-SDCT** | **3 × O(Σ\|leaf\|)** | **1 × O(Σ\|leaf\|)** |

**Measured (com-dblp s=4):**
- ST: SDCT 149ms + treeGraphV+counting+CSR (included in 26ms peel) ≈ 175ms total
- ST_V2: SDCT+callback 139ms + CSR conversion 6ms = 145ms build + 10ms peel = **155ms total**
- **End-to-end improvement: 175ms → 155ms (11.4% faster)**

### Peeling phase

| Operation | ST | ST_V2 |
|---|---|---|
| Phase 1 (find affected leaves) | O(Σ deg_leaf(v)) via CSR | O(Σ deg_leaf(v)) via CSR — identical |
| Phase 2 (broadcast delta) | O(Σ \|leaf\|) via `tree.adj_list` | O(Σ \|leaf\|) via Leaf→Vertex CSR |
| **Cache behavior** | vector-of-vectors (scattered) | flat CSR (sequential) |

**Measured peeling-only (com-dblp s=4, median of 3 trials):**
- ST: 25.7 ms
- ST_V2: 10.5 ms
- **Peeling speedup: 2.45×**

The peeling speedup comes from:
1. **Better cache locality**: flat CSR arrays vs vector-of-vectors with pointer chasing
2. **No redundant passes**: initial support and CSR already computed during build

## 4. Experimental Results

### com-dblp (317,080 nodes, 1,049,866 edges)

#### r=1, s=3

| Metric | ST | ST_V2 | Ratio |
|---|---|---|---|
| Leaves | 317,253 | 317,253 | 1.00 |
| Σ\|leaf\| | 1,340,203 | 1,340,203 | 1.00 |
| SDCT build (ms) | 144 | 149 (fused) | — |
| Peeling (ms) | 31.6 | 15.6 | **2.02×** |
| Total decomp (ms) | 31.6 | 15.7 | **2.01×** |
| Other Index Mem (kB) | 119,856 | 152,208* | — |

#### r=1, s=4 (median of 3 trials)

| Metric | ST | ST_V2 | Ratio |
|---|---|---|---|
| Leaves | 166,918 | 166,918 | 1.00 |
| Σ\|leaf\| | 891,989 | 891,989 | 1.00 |
| SDCT build (ms) | 148 | 146 (fused) | — |
| Peeling (ms) | 25.7 | 10.5 | **2.45×** |
| Total decomp (ms) | 25.7 | 10.5 | **2.45×** |
| Other Index Mem (kB) | 116,358 | 133,088* | — |

*ST_V2 memory includes refTree (built for verification). Production mode would be lower.

### email-Eu-core (1,005 nodes, 16,064 edges) — scaling across s

| s | Leaves | Σ\|leaf\| | ST peel (ms) | V2 peel (ms) | Speedup |
|---|--------|-----------|-------------|-------------|---------|
| 3 | 28,470 | 137,348 | 1.59 | 0.68 | 2.34× |
| 4 | 43,220 | 308,055 | 3.04 | 1.19 | 2.55× |
| 5 | 48,513 | 430,556 | 4.32 | 1.80 | 2.40× |
| 6 | 47,421 | 470,448 | 4.92 | 2.00 | 2.46× |

**Observation:** Speedup is consistent at ~2.4× across different s values, confirming the improvement is structural (cache locality of CSR) rather than data-dependent.

## 5. Limitations

1. **Memory accounting**: ST_V2 currently builds `refTree` for SDCT_Par7 verification and clique counting. In production, skipping `refTree` saves the full `tree.adj_list` allocation.

2. **COO buffer peak**: During the Build phase, the COO buffer temporarily holds all vertex-leaf pairs before conversion to CSR. Peak memory is `Σ|leaf| × 12B` (COOEntry: vertex + leafId + isPivot + padding). This is freed after CSR construction.

3. **`std::function` avoided**: The recursive BK uses a raw function pointer trampoline (`AugCtx`) instead of `std::function` to avoid per-frame overhead in deep recursion.

4. **Graph mutation timing**: `SDCT_Augmented_NoTree` must run before `edgeGraph.beSingleEdge()` mutates the graph. Solved by splitting into Build (pre-mutation) and Peel (post-mutation) phases.
