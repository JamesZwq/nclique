# (r,s)-Nucleus Decomposition: Optimization Summary

本文档总结本代码库中所有算法变体相对于各自 SOTA baseline 的理论与工程优化。

**SOTA Baselines:**
- **R=1**: `NCliqueVertexCoreDecomposition` — 可变树 + Boost d-ary heap 逐顶点剥离
- **R=2**: `PlusNucleusEdgeCoreDecompositionSet` — 可变树 + OMP 并行 + 逐边剥离
- **R≥3**: `NucleusCoreDecompositionCorrect` — 可变树 + BK 重枚举 + 逐 r-clique 剥离

---

## 一、核心基础设施优化：SDCT_Augmented

**文件**: `src/SDCT_Augmented.h`, `src/SDCT_Augmented.inl`

| 维度 | 说明 |
|------|------|
| **核心思想** | 回调式 SDCT：在 BK 递归中每发现一个叶节点，立刻调用 `onLeaf(leafId, keepV, dropV)` 回调，将下游处理（支持度计数、CSR 构建等）融合进单次 SDCT 遍历 |
| **理论优化** | 消除传统两遍工作流（先枚举全部极大团、再后处理）中的冗余遍历；将 3 个 O(Σ\|leaf\|) 后处理 pass 融合为 1 个 |
| **工程优化** | `AugCtx` 函数指针蹦床（trampoline）实现零开销类型擦除——每层递归栈帧仅需 16B（函数指针+void*），避免 `std::function` 在深递归中的栈溢出 |
| **适用范围** | 支撑 R=1 ST_V2（无树）、Interleaved 等变体 |

**两种模式:**
- `bkRecurse_WithTree()`: 构建树 + 触发回调（用于需要树的变体）
- `bkRecurse_NoTree()`: 不建树、仅计数 + 回调（用于无树变体，省去 O(T·k) 树存储）

---

## 二、R=1 顶点核分解优化

### Baseline: NCliqueVertexCoreDecomposition

- 可变树 `DynamicGraph<TreeGraphNode>`，每次剥离都物理修改树
- 顶点→叶索引 `DynamicGraphSet<TreeGraphNode>`（robin_hood hash set）
- Boost d-ary heap（arity=8）管理剥离优先级
- **瓶颈**: 树变异开销大（指针追逐、hash set erase）、缓存不友好

### ST: 不可变树 + 整数计数器 + nCr 差分

**文件**: `NCliqueVertexCoreDecompositionST.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 彻底消除树变异：用 `leafRemainPivots[L]`、`leafNeedPivot[L]`、`leafAlive[L]` 三个整数计数器替代可变树状态。支持度变化量通过 nCr 差分公式 `nCr[old_rp][np] - nCr[new_rp][np]` O(1) 计算，无需 BK 重枚举 |
| **工程优化** | 用扁平 CSR（`vtxLeafOff[]` + `vtxLeafData[]`）替代 hash set 索引；桶数组（bucket array）替代 Boost heap（O(1) 插入/弹出）；CSR 区域预取（prefetch） |
| **复杂度** | 剥离阶段 O(V + E')，E' = CSR 总条目数。Baseline 为 O(V·log V + E'·log V) |
| **加速比** | ~1.5x–2.0x vs Baseline |

### ST_V2: 无树（Tree-Free）

**文件**: `NCliqueVertexCoreDecompositionST_V2.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 完全消除树存储：通过 `SDCT_Augmented_NoTree` 回调式枚举，在 SDCT 遍历中在线构建 dual CSR，从不分配 `tree.adj_list`。省去 ~53B/条目 → ~16B/条目（3.3x 内存削减） |
| **工程优化** | 分离 Build/Peel 两阶段（Build 在 `beSingleEdge()` 之前运行，需要原始图；Peel 在之后运行）；COO → CSR 排列在回调中在线完成 |
| **复杂度** | 与 ST 相同 O(V + E')，但常数因子更低（无树开销） |
| **为何仅 R=1** | R=1 剥离永远不需要 BK 重枚举（只更新计数器），不可变 CSR 即可满足。R≥2 有 Case B 需要动态树变异 |
| **加速比** | 2.3x–4.8x vs ST（视图规模而定） |

### Local (V1–V4): H-index 迭代收敛

**文件**: `NCliqueVertexCoreDecompositionLocal.cpp` (V1), `...V2.cpp`, `...V3.cpp`, `...V4.cpp`

| 变体 | 核心优化 | 并行 |
|------|---------|------|
| **Local V1** | 每个顶点独立计算 H-index（基于相邻叶的支持度），脏队列驱动收敛。整数桶 H-index O(k) 替代 O(k log k) 排序；合并单次叶扫描 | 否 |
| **Local V2** | + core-level 入队过滤（v 的旧 core < u 的 core 则跳过）+ 叶时间戳跳过（叶未变则不重算） | 否 |
| **Local V3** | 同步轮次并行（OpenMP）：Phase 1 并行 H-index 计算，Phase 2 顺序更新 + 收集脏邻居。双缓冲 `newCoreV[]` 避免竞争 | 是（同步） |
| **Local V4** | 异步原地更新：利用 `coreV[]` 单调递减性质，线程直接写入无需同步。每线程本地工作列表，避免原子队列 | 是（异步） |

**理论基础**: H-index 收敛保证正确性（单调递减 → 不动点 = 核值）。V4 利用单调性保证异步安全——线程读到旧值（更高）不影响正确性，最多延迟一轮收敛。

### Interleaved: 构建-分解交错

**文件**: `NCliqueVertexCoreDecompositionInterleaved.cpp`, `...V2.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 在 SDCT 构建过程中同步剥离：当顶点 v 的子树完成时，若有效支持度 ≤ 当前最小核值，立即剥离 v。双层支持度追踪：`countingV[u]`（已完成叶累积）+ `pendingDelta[u]`（已剥离邻居扣减） |
| **工程优化** | 向量式邻接（`vtxLeafAdj[]`, `leafVtxAdj[]`）在回调中增量构建；动态桶数组 |
| **优势** | 若大量顶点在构建期间即可剥离至 0，可减少峰值堆大小 |
| **V2 修复** | 修正 `drainHeap` 逻辑（V1 从未在构建期间真正排空桶）；增加排空守卫；预分配 `vtxLeafAdj` |

### OnDemand: 最小元数据

**文件**: `NCliqueVertexCoreDecompositionOnDemand.cpp`

| 维度 | 说明 |
|------|------|
| **工程优化** | 将 4 个 per-leaf 向量（`leafPivotCount`, `leafNeedPivot`, `leafAlive`, `leafRemainPivots`）精简为 2 个（`leafOrigKeepCount`, `leafOrigPivotCount`），`leafAlive` 状态在需要时即时计算。Init 后释放原始计数向量，节省 ~12B/leaf |

---

## 三、R=2 边核分解优化

### Baseline: PlusNucleusEdgeCoreDecompositionSet (ST 单线程版)

- 可变树 + `treeGraphV`（vertex→leaf hash set）
- **三类情况**:
  - **Case A** (~82–94%): 叶死亡（keepC 被移除）→ 扣除全部贡献
  - **Case B** (~0.5–6.7%): keepC 间边被移除 → BK 重枚举产生子叶
  - **Case C** (~2–10%): 仅 pivot 被移除 → nCr 差分公式
- **Phase 1**: `intersect_dense_sets(treeGraphV[u], treeGraphV[v])` 找受影响叶
- **Phase 2**: O(|leaf|²) 对枚举 × `getEdgeCompressedId` hash 查找 → **主要瓶颈**
- **典型耗时** (facebook s=4): Phase1=328ms, Phase2=1723ms, CaseB=91ms, **Total=2534ms**

### ST_V5: 排序向量替代 hash set

**文件**: `NCliqueEdgeCoreDecompositionPlusSetST_V5.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | Phase 1 交集从 hash probe O(min(\|A\|,\|B\|)) 变为归并扫描 O(\|A\|+\|B\|)，顺序访存 |
| **工程优化** | `SortedVecIndex`：每顶点 `std::vector<TreeGraphNode>`（8B/条目 vs hash set ~60B/条目），脏标记延迟重排序。7.5x 内存削减 |
| **权衡** | 需要 O(n log n) 重排序开销（摊销）；适合内存受限场景 |

### ST_V6: leafAlive 位向量跳过 removeNbr

**文件**: `NCliqueEdgeCoreDecompositionPlusSetST_V6.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | Case A 不再物理删除 hash set 条目：设 `leafAlive[leafId]=0`，在 Phase 1 回调中 `if (!leafAlive[leafId]) return` 过滤。省去 removeNbr 的 hash erase 开销 |
| **失败原因** | 残留条目（stale entries）使 Phase 1 交集迭代大量死条目。facebook s=4: Phase1 从 338ms 恶化至 3552ms（10.5x），Phase2 仅省 176ms。**总体 2x 慢** |

### ST_V7/V7b/V7c: 三种残留条目清理策略

| 变体 | 清理策略 | Phase 1 处理 | 结果 |
|------|---------|-------------|------|
| **V7** (Lazy Purge) | `staleCount[v]` 按需清理：交集前若 `staleCount[v]>0` 则清扫 | 无 leafAlive 检查（已清洁） | Phase1 改善但 Phase2 无变化，总体未胜出 |
| **V7b** (Periodic) | `dirtyVertices` + 每 100 轮批量清扫 | 保留 leafAlive 检查（两次清扫间可能有残留） | Phase1 最佳但清扫本身耗时，总体最慢 |
| **V7c** (Always-Purge) | 每次交集前无条件清扫双方 set | 无 leafAlive 检查 | 最简单但清扫过度，仍未胜出 |

**结论**: leafAlive 路线对 R=2 根本不可行——省下的 Phase2 远不抵 Phase1 退化。

### ST_V8: Case A LeafEdgeInfo 快速路径

**文件**: `NCliqueEdgeCoreDecompositionPlusSetST_V8.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | Case A 从 O(\|leaf\|²) 对枚举 + hash 查找 → O(\|leafEdgeInfo[li]\|) 线性扫描。Init 时预建 `leafEdgeInfo[leafId] = [{edgeId, EdgeType}]` + per-leaf 权重 `leafWKK/PP/KP`。Case A 仅需按类型乘权重 |
| **工程优化** | 一次 init 同时建 `countingKE` 和 `leafEdgeInfo`；`switch(entry.type)` 分派权重；Case C/B 使失效（clear+shrink_to_fit）；Overflow 叶（Case B 子叶）回退原始枚举 |
| **瓶颈** | `vector<vector<LeafEdgeEntry>>` 的 46.6M 条目 init 分配耗时 +515ms |
| **加速比** | Phase2: 1723→1480ms (14%)，但 init 开销大，总体仅 1.04x |
| **Case A 命中率** | 95.2%（1,520,707 fast / 76,652 fallback） |

### ST_V8b: CSR Init + deltaAccum 批量 bucketMove (**当前最优**)

**文件**: `NCliqueEdgeCoreDecompositionPlusSetST_V8b.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化 1 (CSR Init)** | 用扁平 CSR (`leafEdgeOff[numLeaves+1]` + `leafEdgeData[totalEntries]`) 替代 `vector<vector<>>`。两遍构建：Pass 1 计数、前缀和；Pass 2 填充。单次分配替代数百万小向量分配 |
| **理论优化 2 (deltaAccum)** | `directSub` 不再立即调用 `bucketMove`：仅标记 `edgeDirty[id]=1` 并加入 `dirtyEdges`。Phase 2 全部完成后一次性 flush——消除同一边被多个叶多次 bucketMove 的冗余 |
| **工程优化** | `leafEdgeAlive[li]` 追踪 CSR 数据有效性；flush 阶段独立计时 |

**性能分解** (facebook_combined, s=4, 3 次中位数):

| 指标 | ST baseline | V8 | **V8b** |
|------|-----------|-----|---------|
| **Total** | **2534 ms** | 2443 ms | **2048 ms** |
| Init | ~50 ms | 515 ms | 405 ms |
| Phase 1 | 328 ms | 334 ms | 331 ms |
| Phase 2 (A+C) | 1723 ms | 1480 ms | **1134 ms** |
| Case B | 91 ms | 92 ms | 64 ms |
| Flush | — | — | 90 ms |

**CSR Init 改善**: 515→405ms（-21%，消除百万次小向量分配）
**deltaAccum 改善**: Phase2 1480→1134ms（-23%，批量去重 bucketMove）
**总加速比**: **1.24x** vs ST baseline

### Lazy: 按需重算

**文件**: `NCliqueEdgeCoreDecompositionPlusSetST_Lazy.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 不计算精确 Δ：Case A/C 仅标记受影响边为 dirty。边从桶弹出时才按需重算支持度（通过 treeGraphV 交集）。适合大部分边早期被移除的场景 |
| **权衡** | 若边存活至高桶，重算开销 > 预计算 Δ 开销 |

### TreeFree: 消除可变树

**文件**: `NCliqueEdgeCoreDecompositionTreeFree.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 完全消除 `tree.adj_list` 变异：用不可变 dual CSR（`vtxLeafOff/Data` + `leafVtxOff/Data`）替代。Phase 1 用 CSR 归并扫描替代 hash set 交集 |
| **工程优化** | `leafAlive[]` + `leafMutated[]` 追踪状态；Case B 回退到从 `leafVtxData` 重建叶顶点列表 + BK |

### TreeFreeV2: Edge→Leaf 反向索引

**文件**: `NCliqueEdgeCoreDecompositionTreeFreeV2.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 移除 TreeFree 中未使用的 dual CSR（-30% init 内存）；新增 Edge→Leaf CSR (`edgeLeafOff/Data`)，Phase 1 直接从边查叶 O(leaves_per_edge)，无需双端顶点交集 |
| **工程优化** | Case C 后可重建 leafEdgeInfo（小开销），使后续 Case A 仍可走快速路径；`leafMutated[]` 区分 CSR 有效/失效 |

### OnDemand: 最小改动

**文件**: `NCliqueEdgeCoreDecompositionOnDemand.cpp`

| 维度 | 说明 |
|------|------|
| **定位** | ST 的保守变体，Phase 1/2 逻辑基本不变，用于调试对照 |

---

## 四、R≥3 通用 r-Clique 核分解优化

### Baseline: NucleusCoreDecompositionCorrect

- 可变树 + `enumerateCombinations()` + `byClique()` hash map 查找
- 每次叶更新需重新枚举 r-clique、hash 查找所属 clique ID → 高开销
- BK 回调中 hash map 包含性检查 ~30-60 cycles/次

### ST: LeafDeath 快速路径 + BK 分类

**文件**: `NucleusCoreDecompositionRCliqueST.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | `vertexConflictDeg` 分析识别叶死亡（跳过 BK）；仅对"真冲突"叶调用 BK |
| **工程优化** | 桶数组替代 heap；整数算术替代浮点 |

### ST_V2–V4: 预建索引

| 变体 | 核心优化 | 瓶颈 |
|------|---------|------|
| **V2** | `leafCliqueInfo[leafId] = [{cliqueId, ncrValue}]`（叶→clique 反向索引） | 消除 Support 阶段重复 `enumerateCombinations` + hash 查找 |
| **V3** | `cliqueLeafIds[cliqueId] = [leafId]`（clique→叶 反向索引） | 消除 Intersect 阶段 r-way hash set 多路交集 |
| **V4** | V2+V3 融合：单次 init 同时建双向索引 | 更好的缓存局部性 + 减少遍历次数 |

### ST_V5–V8: BK 加速 + Case C 提取

| 变体 | 核心优化 | 提升方向 |
|------|---------|---------|
| **V5** | 位置化包含检查：`c.test(pos)` 替代 hash 查找（1 cycle vs 30-60 cycles） | BK 回调 20-40x 加速 |
| **V6** | Case C 提取：纯 pivot 冲突用 nCr 差分公式，跳过 BK | BK 调用率降低 85-97% |
| **V7** | 放宽 Case C：所有冲突顶点均为 pivot 即可（不要求满冲突度） | 进一步降低 BK 调用率至 <5% |
| **V8** | BK 回调用 bitset + `vListMap` 位置映射，彻底消除 hash map | BK 包含检查 30-50x 加速 |

### ST_V9–V11: 内存与微优化

| 变体 | 核心优化 | 提升方向 |
|------|---------|---------|
| **V9** | 消除剥离阶段 treeGraphV 维护（只写不读）；`coverToVertex` 缓冲区复用 | 内存 -30-50%，减少 malloc |
| **V10** | BK 回调微优化：`newEntries` 复用（避免 per-callback 分配）；`addNodePresorted()`（跳过已排序的 sort） | BK 回调 ~5-10% 加速 |
| **V11** | `uint64_t posMask` 位掩码：`(posMask & childMask) == posMask`（1 AND + 1 CMP，无分支）；`popcnt` 计算 pivot 数 | 24B/条目 vs 32B（-25%内存）+ 极快包含检查 |

### ST_V12: 按需重枚举（内存优化）

**文件**: `NucleusCoreDecompositionRCliqueST_V12.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 完全消除 `leafCliqueInfo`：剥离时按需重枚举 r-clique（`enumerateCombinationsWithIdx` + `lookupRaw`）。per-leaf 缓存在叶处理后释放 |
| **空间** | O(\|C\| + \|V\|) vs baseline O(K × C(avg_leaf, r))。**节省 60-80% 内存** |
| **权衡** | 每叶重枚举增加 CPU 开销（~10-15%），但对大图内存受限场景至关重要 |

### ST_V13: 路径脆弱性剥离（PFP）

**文件**: `NucleusCoreDecompositionRCliqueST_V13.cpp`

| 维度 | 说明 |
|------|------|
| **理论优化** | 彻底消除 `StaticCliqueIndex` 和所有 r-clique 枚举。按"路径脆弱性"（path fragility）剥离：孤立路径（~85%+）用 O(r) 闭式公式计算核值；非孤立路径按需顶点-路径交集 |
| **空间** | O(V + Σ\|P\|) — 仅 CPI 本身大小，无 r-clique 枚举开销 |
| **适用场景** | 高孤立度工作负载（稀疏或聚类图） |

---

## 五、优化分类总览

### 理论优化（算法层面改进）

| 优化 | 适用 | 描述 |
|------|------|------|
| **不可变树 + nCr 差分** | R=1 ST, R=2 ST | 消除树变异，支持度变化量 O(1) 闭式计算 |
| **无树架构** | R=1 ST_V2 | 回调式 SDCT 在线构建 CSR，省去 O(T·k) 树存储 |
| **H-index 迭代收敛** | R=1 Local | 异步局部计算，天然适合并行化 |
| **Case 分类** | R=2 ST, R≥3 V6/V7 | Case A/C 用闭式公式，仅 Case B 需 BK |
| **LeafEdgeInfo 索引** | R=2 V8/V8b | Case A 从 O(\|leaf\|²) 降至 O(\|leafEdgeInfo\|) |
| **deltaAccum 批量 bucketMove** | R=2 V8b | 消除同一边多次冗余桶移动 |
| **按需重枚举** | R≥3 V12 | 用 CPU 换内存：不存 leafCliqueInfo，按需重算 |
| **路径脆弱性** | R≥3 V13 | 消除 r-clique 枚举，用路径级闭式公式 |

### 工程优化（数据结构与实现层面）

| 优化 | 适用 | 描述 |
|------|------|------|
| **扁平 CSR 替代 hash set** | R=1 ST/V2, R=2 TreeFree/V2 | 顺序访存、单次分配、缓存友好 |
| **桶数组替代 heap** | 全部 ST 变体 | O(1) 弹出/插入 vs O(log n) heap |
| **位掩码包含检查** | R≥3 V5/V8/V11 | 1 AND + 1 CMP 替代 hash 查找（30-50x） |
| **函数指针蹦床** | SDCT_Augmented | 深递归中避免 std::function 栈溢出 |
| **两遍 CSR 构建** | R=2 V8b | Pass1 计数 + 前缀和，Pass2 填充，单次分配 |
| **排序向量归并交集** | R=2 V5/TreeFree | O(\|A\|+\|B\|) 顺序扫描替代 hash probe |
| **Edge→Leaf 反向 CSR** | R=2 TreeFreeV2 | Phase1 直接从边查叶，无需双端交集 |
| **per-thread 暂存缓冲** | R=1 LocalV3/V4 | 消除并行竞争，每线程独立 scratch buffer |
| **单调性保证异步安全** | R=1 LocalV4 | `coreV[]` 只降不升 → 线程读到旧值仍正确 |

---

## 六、性能实测汇总

### R=1: ST_V2 vs ST (单线程)

| 图 | 顶点数 | s=3 | s=4 | s=5 |
|----|-------|-----|-----|-----|
| com-dblp | 317K | 2.3x | 2.7x | 1.5x |
| com-youtube | 1.13M | 3.4x | 4.1x | 4.8x |
| web-Google | 876K | 2.4x | 2.4x | 2.7x |
| soc-pokec | 1.63M | 2.3x | 2.7x | 3.4x |

### R=2: V8b vs ST baseline (单线程, 3 次中位数)

| 图 | s | ST | V8b | 加速比 |
|----|---|-----|------|--------|
| facebook_combined | 3 | 236 ms | 216 ms | 1.09x |
| facebook_combined | 4 | 2534 ms | 2048 ms | **1.24x** |
| com-dblp | 3 | 384 ms | 359 ms | 1.07x |
| com-dblp | 4 | 251 ms | 234 ms | 1.07x |

### R≥3: ST_V11 vs Correct

| 图 | 加速比范围 |
|----|----------|
| 多图综合 | 1.7x–3.2x |

---

## 七、各 r 值当前最优变体

| r | 最优变体 | 核心优势 |
|---|---------|---------|
| **R=1** | **ST_V2** | 无树 + 回调融合 CSR，2.3–4.8x 加速 |
| **R=2** | **V8b** | CSR init + deltaAccum，1.24x 加速（Phase2 主导场景） |
| **R≥3** | **ST_V11** | 位掩码 + 双向索引 + Case C 提取，1.7–3.2x 加速 |
| **R≥3 (内存受限)** | **ST_V12** | 按需重枚举，60-80% 内存节省 |
