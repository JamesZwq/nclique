# 冲突分析：tree/treeGraphV 更新的真实情况

## 场景分析

假设我们同时处理两个 leaves：L1 和 L2

### Leaf L1 的更新
1. 移除 L1 from tree
2. 对于 L1 中的每个顶点 v：从 treeGraphV[v] 中移除 L1
3. 添加新 leaves L1_new1, L1_new2, ...
4. 对于新 leaves 中的每个顶点 v：向 treeGraphV[v] 添加新 leaf

### Leaf L2 的更新
1. 移除 L2 from tree
2. 对于 L2 中的每个顶点 v：从 treeGraphV[v] 中移除 L2
3. 添加新 leaves L2_new1, L2_new2, ...
4. 对于新 leaves 中的每个顶点 v：向 treeGraphV[v] 添加新 leaf

## 冲突点分析

### 冲突点 1：tree 的全局 ID 分配
- `tree.addNode()` 需要分配新的 leaf ID
- 这是一个全局计数器
- **解决方案**：使用 atomic counter

### 冲突点 2：treeGraphV[v] 的并发修改
- 如果 L1 和 L2 都包含顶点 v
- 两个线程会同时修改 treeGraphV[v]
- **这是真正的冲突！**

但是，等等...

## 关键洞察：操作的原子性

### 当前的操作
```cpp
treeGraphV.removeNbr(v, {leafId, isPivot});  // 从 vector 中删除
treeGraphV.addNbr(v, {newLeafId, isPivot});  // 向 vector 中添加
```

### 问题
- `removeNbr` 和 `addNbr` 都修改 `treeGraphV[v]` 这个 vector
- 如果两个线程同时操作，vector 会损坏

### 解决方案 1：细粒度锁（每个顶点一个锁）
```cpp
std::vector<omp_lock_t> vertexLocks(graphN);

// 初始化
for (auto &lock : vertexLocks) omp_init_lock(&lock);

// 使用
omp_set_lock(&vertexLocks[v]);
treeGraphV.removeNbr(v, {leafId, isPivot});
omp_unset_lock(&vertexLocks[v]);
```

**问题**：锁的开销很大，而且可能死锁

### 解决方案 2：无锁队列
```cpp
// 每个顶点维护两个队列
std::vector<std::vector<std::pair<LeafId, bool>>> toRemove(graphN);
std::vector<std::vector<std::pair<LeafId, bool>>> toAdd(graphN);

// 并行阶段：收集操作
#pragma omp parallel for
for (each leaf) {
    for (v in leaf) {
        toRemove[v].push_back({leafId, isPivot});
    }
    for (v in newLeaf) {
        toAdd[v].push_back({newLeafId, isPivot});
    }
}

// 串行阶段：批量应用
for (v = 0; v < graphN; ++v) {
    for (auto op : toRemove[v]) {
        treeGraphV.removeNbr(v, op);
    }
    for (auto op : toAdd[v]) {
        treeGraphV.addNbr(v, op);
    }
}
```

**这个方案可行！**

### 解决方案 3：使用并发数据结构
```cpp
// 使用 concurrent_vector 或类似的数据结构
tbb::concurrent_vector<TreeGraphNode> treeGraphV[graphN];

// 可以安全地并行添加和删除
#pragma omp parallel for
for (each leaf) {
    for (v in leaf) {
        // 并发安全
        treeGraphV[v].remove({leafId, isPivot});
        treeGraphV[v].add({newLeafId, isPivot});
    }
}
```

## 最有希望的方案：解决方案 2（批量操作）

### 优点
1. 完全并行收集操作
2. 串行应用很快（只是简单的 vector 操作）
3. 无锁，无死锁风险
4. 实现简单

### 预期性能
- 并行收集：O(changedLeaf.size() * leafSize) / nthreads
- 串行应用：O(graphN * avgOperations)

如果 changedLeaf.size() 很大，并行收集的收益会很显著！

## 让我实现这个方案！
