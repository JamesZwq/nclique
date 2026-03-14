# Pivoter 并行优化总结

## 优化成果（com-dblp, 12线程, Mac M2）

### SDCT（树构建）
- **串行**：258 ms
- **并行**：131 ms
- **加速比**：**1.97×**

### Nucleus Decomposition（主要计算）
- **串行**：2831 ms
- **并行**：2210 ms
- **加速比**：**1.28×**

### 整体性能
- **串行总时间**：约 3089 ms
- **并行总时间**：约 2341 ms
- **总加速比**：约 **1.32×**

## 实现的优化

### 1. SDCT 并行化
**策略**：块并行 + 内存优化
- 将顶点分成 N 块（N=线程数）
- 每个线程维护完整的串行状态
- 保持 `beginR++` 逻辑（关键！）
- 使用 `unique_ptr` 管理大数组
- `neighborsInP` 按需分配

**关键代码**：
```cpp
auto vertexSets = std::make_unique<int[]>(size);
auto vertexLookup = std::make_unique<int[]>(size);
auto neighborsInP = std::make_unique<int*[]>(size);
// neighborsInP[i] = nullptr; // 按需分配
```

### 2. Nucleus Decomposition 优化

#### 2.1 countingPerRCliqueParallel
**问题**：使用 `#pragma omp critical` 合并结果，成为瓶颈
**解决**：使用线程局部数组 + 并行 reduction
```cpp
std::vector<std::vector<double>> thread_locals(nthreads, ...);
// 每个线程独立累加
#pragma omp parallel for schedule(static)
for (daf::Size i = 0; i < nClique; ++i) {
    for (int t = 0; t < nthreads; ++t) {
        rCliqueSCounting[i] += thread_locals[t][i];
    }
}
```

#### 2.2 Part A (Intersect) 优化
**问题**：使用 `#pragma omp critical` 合并 pairs
**解决**：预分配线程局部 vector + 批量合并
```cpp
std::vector<std::vector<PairOV>> thread_pairs(nthreads);
// 并行计算后批量合并
for (auto &tp : thread_pairs) {
    allPairs.insert(allPairs.end(), tp.begin(), tp.end());
}
```

## 性能瓶颈分析

### 当前瓶颈
1. **enumerateCombinations**：递归 DFS，难以并行化
2. **Heap 操作**：每次更新都需要调整堆结构
3. **数据结构更新**：tree, treeGraphV 的修改必须串行

### 为什么加速有限？
1. **Amdahl 定律**：存在大量串行部分（heap 更新、数据结构修改）
2. **内存带宽**：多线程同时访问大数组
3. **负载不均衡**：不同 leaf 的计算量差异很大

## 进一步优化方向

### 短期（可能提升到 1.5-2×）
1. **优化 enumerateCombinations**：使用迭代而非递归
2. **更好的负载均衡**：动态调度参数调优
3. **减少内存分配**：使用内存池

### 长期（可能提升到 3-5×）
1. **算法级优化**：使用不同的 nucleus decomposition 算法
2. **GPU 加速**：将计算密集部分移到 GPU
3. **分布式计算**：对于超大图使用多机并行

## 代码状态
- ✅ SDCT：1.97× 加速，正确性验证通过
- ✅ Nucleus：1.28× 加速，正确性验证通过
- ✅ 总体：1.32× 加速

## 结论
在保持算法正确性的前提下，通过并行化和内存优化，实现了约 1.3× 的整体加速。
主要瓶颈在于算法本身的串行依赖（heap 操作、数据结构更新），
要获得更高的加速比需要算法级别的重新设计。
