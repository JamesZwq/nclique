# Nucleus Decomposition: All Versions

Baseline for all comparisons: **Correct** (`NucleusCoreDecompositionCorrect`), which uses a d-ary heap, full BK enumeration, and `double` support values.

---

## r=1 Versions (Vertex Core Decomposition)

### Correct (Baseline)
- **File**: `NUcllearCoreDeompositoinCorrect.cpp`
- **Env**: `PIVOTER_RUN_REF=1`
- **Data structure**: d-ary heap (boost), `double` support, full BK + r-clique hash index
- **Complexity**: Each peeling step does BK enumeration on affected leaves, enumerates all 1-cliques (vertices), looks up in hash index, updates heap

### ST (Single-Thread Peeling)
- **File**: `NCliqueVertexCoreDecompositionST.cpp`
- **Env**: `PIVOTER_RUN_ST=1`
- **Innovations vs Correct**:
  1. **No r-clique hash index**: r=1 means r-cliques are individual vertices — directly index by vertex ID, eliminating the entire `StaticCliqueIndex` hash table
  2. **Integer counters**: `long long` instead of `double`, avoids floating-point rounding
  3. **Bucket array**: O(1) insert/decrease-key replaces O(log n) d-ary heap
  4. **No BK enumeration**: When a vertex is peeled from a leaf, the leaf simply loses that vertex — no BK needed. Support delta computed analytically via nCr
  5. **Direct tree modification**: Remove vertex from leaf → update nCr contribution to remaining vertices. No new sub-leaf creation
- **Speedup**: ~5-10x vs Correct on large graphs

### Local V1 (H-index Iterative Convergence)
- **File**: `NCliqueVertexCoreDecompositionLocal.cpp`
- **Env**: `PIVOTER_RUN_LOCAL=1`
- **Innovations vs Correct**:
  1. **No peeling / no tree modification**: Core values computed via iterative H-index convergence. The SDCT tree is never modified — read-only after construction
  2. **H-index definition**: `core(v) = max c such that v participates in ≥ c s-cliques where all other vertices have core ≥ c`
  3. **Merged leaf-scan pass**: Single pass over each leaf's nodes to compute minKeepCore and collect pivotCores simultaneously (vs two separate passes)
  4. **Integer bucket H-index**: Bucket array for H-index computation — O(B + maxLevel) instead of O(B log B) sort on breakpoints
  5. **Early termination**: Pass current `coreV[v]` into H-index computation; once accumulated support ≥ currentCore, return immediately
  6. **Dirty queue**: Only re-evaluate vertices whose neighbors' cores decreased (FIFO queue with dedup)
- **Trade-off**: No tree modification means lower memory overhead and simpler code, but may require more iterations to converge

### Local V2 (Core-Level Enqueue Filter + Timestamp Skip)
- **File**: `NCliqueVertexCoreDecompositionLocal_V2.cpp`
- **Env**: `PIVOTER_RUN_LOCAL_V2=1`
- **Innovations vs Local V1**:
  1. **Core-level enqueue filter**: When vertex v's core drops from oldCore, only enqueue neighbor u if `oldCore >= coreV[u]` — vertices with much higher cores can't be affected
  2. **Leaf-level timestamp skip**: Track last-modified timestamp per leaf; skip leaf processing if no participating vertex has changed since last evaluation
- **Speedup**: ~25% faster than V1 on large graphs

### Local V3 (Synchronous Parallel)
- **File**: `NCliqueVertexCoreDecompositionLocal_V3.cpp`
- **Env**: `PIVOTER_RUN_LOCAL_V3=1`
- **Innovations vs Local V2**:
  1. **Round-based parallel H-index**: All queued vertices compute H-index in parallel (reading frozen coreV), then updates applied sequentially
  2. **Per-thread scratch buffers**: pivotCores and buckets arrays are per-thread — no contention
  3. **Parallel Phase 1**: Initial full scan over all vertices parallelized with `omp for`
- **Trade-off**: Synchronous rounds mean within-round updates are invisible → more iterations needed (roughly 2x re-evaluations vs V2)

### Local V4 (Asynchronous Parallel)
- **File**: `NCliqueVertexCoreDecompositionLocal_V4.cpp`
- **Env**: `PIVOTER_RUN_LOCAL_V4=1`
- **Innovations vs Local V3**:
  1. **Async in-place updates**: Threads write `coreV[]` directly (no double buffer). Safe because coreV is monotone decreasing — any stale read is an upper bound
  2. **Per-thread local work lists**: Each thread collects dirty neighbors into its own list, merged after each round. Eliminates sequential neighbor collection bottleneck
  3. **Core-level enqueue filter**: Restored (was removed in V3 due to interaction with sequential apply). Works correctly with async model
  4. **Atomic CAS dedup**: `compare_exchange_strong` on per-vertex `inQueue` flags for O(1) dedup
- **Speedup**: 3.2x–7.7x over V3 sync at 32 threads; 12.4x parallel scaling on soc-pokec (1T→64T)

---

## r=2 Versions (Edge Core Decomposition)

### Correct (Baseline)
- Same as r=1 Correct, but r-cliques are edges → uses edge hash index
- Full BK enumeration for every affected leaf

### ST (Single-Thread Peeling)
- **File**: `NCliqueEdgeCoreDecompositionPlusSetST.cpp`
- **Env**: `PIVOTER_RUN_ST=1`
- **Innovations vs Correct**:
  1. **Case A/B/C classification** for each affected leaf:
     - **Case A (Leaf Death, ~94%)**: Removed edge contains a keep-keep pair, or remaining pivots < needPivot → leaf is destroyed. Full contribution subtracted analytically. **No BK needed.**
     - **Case B (BK Required, ~0.6-5.6%)**: Removed pivot-pivot edges create internal conflicts → must run BK to find surviving sub-leaves
     - **Case C (Pivot-Only Removal, ~1-5%)**: Only pivot vertices removed → analytical delta via `nCr[old] - nCr[new]`. **No BK needed.**
  2. **Integer arithmetic**: `long long countingKE` eliminates float drift
  3. **Bucket array**: O(1) decrease-key
  4. **Immediate bucket move**: No deferred dirty tracking
  5. **Merged Phase 2a+2b**: Delta computation + tree mutation in single pass (avoids re-iterating removedLeaf)
  6. **No OMP overhead**: No locks, no per-thread vectors
- **Speedup**: ~2x vs Correct; Case A dominates (94%), effectively eliminating BK for almost all leaves

### Plus Parallel
- **File**: `NCliqueEdgeCoreDecompositionPlusSet.cpp`
- **Env**: (default for r=2)
- **Innovations vs Correct**:
  1. **Same Case A/B/C classification** as ST
  2. **Parallel Phase 1**: Edge→leaf intersect parallelized with `omp parallel for`
  3. **Parallel Phase 2a**: Case A and Case C delta computations parallelized
  4. **Serial Phase 2b**: BK for Case B leaves + tree mutations remain serial (data-dependent)
  5. **Parallel support counting**: Per-thread local arrays merged after counting
- **Speedup**: ~4x vs Correct at 32 threads

---

## r≥3 Versions (General r-Clique Core Decomposition)

### Correct (Baseline)
- **File**: `NUcllearCoreDeompositoinCorrect.cpp`
- **Env**: `PIVOTER_RUN_REF=1`
- **Data structure**: d-ary heap, `double` support, full BK + `StaticCliqueIndex` hash
- For every affected leaf: BK enumeration (`bkRmClique::removeRClique`), enumerate all r-cliques in new sub-leaves, hash lookup to update support

### ST (Single-Thread Optimized)
- **File**: `NucleusCoreDecompositionRCliqueST.cpp`
- **Env**: `PIVOTER_RUN_ST=1`
- **Innovations vs Correct**:
  1. **Leaf-death fast path**: Before calling BK, check if the leaf can possibly survive. Uses BK's own forced-removal criterion: vertex v is forced out if `conflictDeg(v) >= C(n-1, r-1)` (appears in ALL possible r-cliques containing it). If a keep vertex is forced out → leaf dies. If too few pivots remain → leaf dies. **Skips BK entirely for 17-85% of leaves.**
  2. **Integer arithmetic**: `long long` with `(long long)(x + 0.5)` rounding instead of `double`
  3. **Bucket array**: O(1) decrease-key replaces O(log n) d-ary heap
  4. **No OMP overhead**: No thread_pairs, no parallel merge/sort, no atomic operations. Single-thread code path has no synchronization cost.
  5. **Immediate bucket move**: Bucket position updated inline during support subtraction
  6. **Merged phases**: Intersect → BK → tree mutation → support update in single pass per leaf (vs separate Phase A/B/C in parallel version)
  7. **Serial counting**: Dedicated serial `countingPerRClique` avoids parallel overhead for single-thread case
- **Speedup**: 1.1x–1.5x vs Correct (depends on leaf-death rate)

### RClique Parallel (Existing)
- **File**: `NucleusCoreDecompositionRemoveSclique.cpp`
- **Env**: (default for r≥3)
- **Innovations vs Correct**:
  1. **Bucket array**: Same O(1) decrease-key as ST
  2. **Parallel Phase A**: r-clique→leaf intersect parallelized with per-thread pair lists + sort + merge
  3. **Parallel Phase B**: BK enumeration + r-clique support computation parallelized per leaf
  4. **Serial Phase C**: Tree mutations + bucket updates remain serial
  5. **Parallel counting**: Per-thread local arrays merged after initial support computation
- **Note**: At 1 thread, slower than ST due to parallel infrastructure overhead (thread_pairs allocation, sort, merge)

---

## Version Comparison Matrix

| Feature | Correct | r=1 ST | r=1 Local V1-V4 | r=2 ST | r=2 Plus | r≥3 ST | r≥3 Parallel |
|---------|---------|--------|-----------------|--------|----------|--------|-------------|
| Heap type | d-ary | bucket | none (H-index) | bucket | bucket | bucket | bucket |
| Arithmetic | double | long long | double | long long | long long | long long | double |
| BK needed | always | never | never | 0.6-5.6% | 0.6-5.6% | 15-83% | always |
| Tree modified | yes | yes | **no** | yes | yes | yes | yes |
| Parallel | no | no | V3/V4 only | no | yes | no | yes |
| Hash index | yes | **no** | **no** | **no** | **no** | yes | yes |

---

## Environment Variables Reference

| Variable | Description |
|----------|-------------|
| `PIVOTER_RUN_REF=1` | Run Correct (reference) version |
| `PIVOTER_RUN_ST=1` | Run ST (single-thread optimized) version |
| `PIVOTER_RUN_LOCAL=1` | Run Local V1 (r=1 only) |
| `PIVOTER_RUN_LOCAL_NAIVE=1` | Run Local Naive (r=1 only, full scan every iteration) |
| `PIVOTER_RUN_LOCAL_V2=1` | Run Local V2 (r=1 only) |
| `PIVOTER_RUN_LOCAL_V3=1` | Run Local V3 parallel (r=1 only) |
| `PIVOTER_RUN_LOCAL_V4=1` | Run Local V4 async parallel (r=1 only) |
| `PIVOTER_COMPARE=1` | Enable correctness comparison (run both optimized + Correct) |
| `OMP_NUM_THREADS=N` | Set thread count for parallel versions |

## Build & Run

```bash
mkdir -p build && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/bin/degeneracy_cliques <graph.edges> <r> <s> degen
```
