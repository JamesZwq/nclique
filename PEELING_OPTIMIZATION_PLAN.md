# Peeling 循环优化计划 - 系统性尝试所有方法

## 🎯 目标
优化主循环（peeling）部分，这是花费时间最多的地方

## 📊 当前瓶颈分析

### 主循环结构
```cpp
while (!heap.empty()) {
    // 1. 从 heap 中取出最小 core 的 r-cliques
    // 2. Phase A: 找到受影响的 leaves（已优化）
    // 3. Phase B: 并行处理每个 leaf
    //    - 使用 Bron-Kerbosch 生成新 leaves
    //    - 计算 support 增量和减量
    // 4. 并行更新 support（已优化）
    // 5. 更新 tree 和 treeGraphV
    // 6. 更新 heap
}
```

### 时间分布（估计）
- Phase A（收集 leaves）: ~10%
- Phase B（处理 leaves）: ~40%
- Support 更新: ~15%
- Tree/treeGraphV 更新: ~30%
- Heap 更新: ~5%

## 🔬 优化方案列表

### 方案 1: 优化 Bron-Kerbosch 算法
**目标**: Phase B 占 40% 时间，优化它能带来显著提升

**尝试 1.1**: 缓存中间结果
```cpp
// 为每个 leaf 缓存 Bron-Kerbosch 的中间状态
std::unordered_map<LeafId, BKState> bkCache;
```

**尝试 1.2**: 使用更快的位运算
```cpp
// 使用 bitset 代替 vector 进行集合操作
std::bitset<MAX_VERTICES> candidates;
```

**尝试 1.3**: 预计算邻接关系
```cpp
// 预先计算 leaf 内部的邻接矩阵
bool adjMatrix[MAX_LEAF_SIZE][MAX_LEAF_SIZE];
```

### 方案 2: 优化 heap 操作
**目标**: 减少 heap 更新的开销

**尝试 2.1**: 批量 heap 更新
```cpp
// 收集所有需要更新的 cliques，一次性更新
std::vector<std::pair<CliqueId, double>> updates;
// ... 收集 ...
for (auto [id, newValue] : updates) {
    heap.update(heapHandles[id]);
}
```

**尝试 2.2**: 使用更快的 heap 实现
```cpp
// 尝试 Fibonacci heap 或 pairing heap
using Heap = boost::heap::pairing_heap<...>;
```

**尝试 2.3**: 延迟 heap 更新
```cpp
// 只标记需要更新，实际更新延迟到下次使用
std::vector<bool> needUpdate(nClique, false);
```

### 方案 3: 优化 support 计算
**目标**: 减少 nCr 查找和计算开销

**尝试 3.1**: 预计算常用的 nCr 值
```cpp
// 为常见的参数组合预计算
std::unordered_map<std::pair<int,int>, double> nCrCache;
```

**尝试 3.2**: 使用 SIMD 加速
```cpp
// 使用 AVX2/AVX512 并行计算多个 support 值
__m256d values = _mm256_load_pd(nCr[i]);
```

**尝试 3.3**: 减少重复计算
```cpp
// 缓存每个 leaf 的 pivotC 和 keepC
struct LeafInfo {
    int pivotC, keepC;
    // ...
};
```

### 方案 4: 优化内存访问模式
**目标**: 提高缓存命中率

**尝试 4.1**: 数据局部性优化
```cpp
// 将相关数据放在一起
struct CliqueData {
    double support;
    bool inHeap;
    HeapHandle handle;
} __attribute__((aligned(64)));
```

**尝试 4.2**: 预取数据
```cpp
// 提前预取下一次迭代需要的数据
__builtin_prefetch(&tree.adj_list[nextLeafId]);
```

**尝试 4.3**: 使用内存池
```cpp
// 为频繁分配的对象使用内存池
MemoryPool<std::vector<TreeGraphNode>> leafPool;
```

### 方案 5: 减少同步开销
**目标**: 进一步减少线程同步

**尝试 5.1**: 无锁队列
```cpp
// 使用无锁队列收集结果
boost::lockfree::queue<LeafResult> resultQueue;
```

**尝试 5.2**: 分区处理
```cpp
// 将 leaves 按某种规则分区，减少冲突
std::vector<std::vector<LeafId>> partitions;
```

**尝试 5.3**: 流水线并行
```cpp
// 将不同阶段流水线化
#pragma omp parallel sections
{
    #pragma omp section
    { /* Phase A */ }
    #pragma omp section
    { /* Phase B */ }
}
```

### 方案 6: 算法层面优化
**目标**: 从算法角度减少工作量

**尝试 6.1**: 增量更新
```cpp
// 只更新真正改变的部分
if (leaf_changed_significantly) {
    // 完全重新计算
} else {
    // 增量更新
}
```

**尝试 6.2**: 剪枝
```cpp
// 跳过不会影响结果的 leaves
if (leaf.maxPossibleCore < currentMinCore) {
    continue;
}
```

**尝试 6.3**: 批处理
```cpp
// 一次处理多个相似的 r-cliques
std::vector<CliqueId> batch;
for (auto id : currentRemoveRcliqueIds) {
    if (similar(id, batch.back())) {
        batch.push_back(id);
    } else {
        processBatch(batch);
        batch.clear();
    }
}
```

### 方案 7: 并行粒度调整
**目标**: 找到最优的并行粒度

**尝试 7.1**: 动态调整 schedule
```cpp
// 测试不同的 schedule 参数
#pragma omp for schedule(dynamic, 1)   // 当前
#pragma omp for schedule(dynamic, 4)   // 尝试
#pragma omp for schedule(dynamic, 16)  // 尝试
#pragma omp for schedule(guided)       // 尝试
```

**尝试 7.2**: 嵌套并行
```cpp
// 在 leaf 处理内部再并行
#pragma omp parallel for num_threads(2)
for (auto &newLeaf : newLeaves) {
    // 计算 support
}
```

**尝试 7.3**: 任务并行
```cpp
// 使用 task 而不是 for
#pragma omp task
{
    processLeaf(leafId);
}
```

## 📋 实验计划

### 阶段 1: 快速测试（每个 5-10 分钟）
1. 方案 7.1: 调整 schedule 参数
2. 方案 2.3: 延迟 heap 更新
3. 方案 3.3: 减少重复计算
4. 方案 4.2: 数据预取

### 阶段 2: 中等改动（每个 30-60 分钟）
1. 方案 1.2: 使用位运算
2. 方案 2.1: 批量 heap 更新
3. 方案 5.2: 分区处理
4. 方案 6.2: 剪枝优化

### 阶段 3: 大改动（每个 2-4 小时）
1. 方案 1.1: 缓存 BK 状态
2. 方案 2.2: 更换 heap 实现
3. 方案 4.3: 内存池
4. 方案 5.1: 无锁队列

## 🔄 实验流程

对于每个方案：
1. 创建新分支/备份
2. 实现优化
3. 编译
4. 在 com-dblp 上快速测试
5. 如果有改进，在 web-Google 上完整测试
6. 记录结果
7. 如果改进 > 5%，保留；否则回滚

## 📊 性能记录模板

```
方案: [方案编号]
描述: [简短描述]
实现时间: [分钟]
编译: [成功/失败]

com-dblp 结果:
- 之前: XXX ms
- 之后: XXX ms
- 改进: XX%

web-Google 结果:
- 之前: XXX ms
- 之后: XXX ms
- 改进: XX%

结论: [保留/回滚]
备注: [其他观察]
```

## 🎯 目标

通过系统性尝试所有方案，找到能够带来显著提升的优化组合。

**预期**：
- 如果每个方案平均提升 5%
- 10 个成功的方案 = (1.05)^10 = 1.63× 额外加速
- 结合当前的 1.69× = 1.69 × 1.63 = **2.75× 总体加速**

**乐观**：
- 如果找到几个关键优化（如 BK 算法优化 20%）
- 可能达到 **3-4× 的 Nucleus 加速**
- 结合 SDCT_Par，总体可能达到 **6-8× 加速**

让我们开始系统性地尝试每一个方案！
