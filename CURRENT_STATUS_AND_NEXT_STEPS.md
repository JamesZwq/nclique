# 当前状态与下一步优化计划

## 当前已实现的优化

### 1. ✅ countingPerRClique 优化
- 移除 `#pragma omp critical`
- 使用线程局部数组 + 并行 reduction
- **效果**：1.68-2.72× 加速

### 2. ✅ Part A 优化（收集 leaves）
- 移除 `#pragma omp critical`
- 使用线程局部缓冲区
- **效果**：有一定改善

### 3. ✅ 并行 support 更新
- 使用 `std::atomic<double>` + `compare_exchange_weak`
- 完全并行更新所有 cliques 的 support
- **效果**：web-Google 显著改善

### 当前最佳性能
- **web-Google（64线程）**：830 ms → 492 ms（**1.69× 加速**）
- **com-dblp（8线程）**：498 ms → 436 ms（**1.14× 加速**）

## 发现的问题

### BatchParallel 算法崩溃
**症状**：程序在处理第一批 r-cliques 后崩溃（Segmentation fault）

**可能原因**：
1. 内存访问越界
2. 数据结构损坏
3. 并发访问冲突

## 下一步优化计划

### 方案 1：修复并完善批量操作（推荐）

#### 步骤 1：简化实现，逐步测试
```cpp
// 最简单的版本：只批量收集，串行应用
std::vector<std::vector<TreeGraphNode>> toRemove(graphN);
std::vector<std::vector<TreeGraphNode>> toAdd(graphN);

// 串行收集（先确保逻辑正确）
for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
    auto leafId = changedLeaf[idx];
    const auto leaf = tree.adj_list[leafId];
    
    for (auto leafV : leaf) {
        toRemove[leafV.v].push_back({leafId, leafV.isPivot});
    }
    
    // ... 更新 tree ...
    
    for (auto &newLeaf : res.newLeaves) {
        daf::Size newId = tree.addNode(newLeaf);
        for (auto i : tree.adj_list[newId]) {
            toAdd[i.v].push_back({newId, i.isPivot});
        }
    }
}

// 批量应用 treeGraphV 更新
for (daf::Size v = 0; v < graphN; ++v) {
    for (const auto &node : toRemove[v]) {
        treeGraphV.removeNbr(v, node);
    }
    for (const auto &node : toAdd[v]) {
        treeGraphV.addNbr(v, node);
    }
}
```

#### 步骤 2：验证正确性
- 运行小数据集
- 对比结果与原版本
- 确保没有崩溃

#### 步骤 3：逐步并行化
```cpp
// 并行批量应用
#pragma omp parallel for schedule(dynamic, 1024)
for (daf::Size v = 0; v < graphN; ++v) {
    // 批量更新 treeGraphV[v]
}
```

#### 步骤 4：测试性能
- 测试不同线程数
- 对比优化前后
- 分析瓶颈

### 方案 2：使用细粒度锁（备选）

```cpp
// 为每个顶点创建一个锁
std::vector<omp_lock_t> vertexLocks(graphN);
for (auto &lock : vertexLocks) omp_init_lock(&lock);

// 并行更新，使用细粒度锁
#pragma omp parallel for schedule(dynamic, 16)
for (daf::Size idx = 0; idx < changedLeaf.size(); ++idx) {
    auto leafId = changedLeaf[idx];
    const auto leaf = tree.adj_list[leafId];
    
    // 更新 treeGraphV（使用锁）
    for (auto leafV : leaf) {
        omp_set_lock(&vertexLocks[leafV.v]);
        treeGraphV.removeNbr(leafV.v, {leafId, leafV.isPivot});
        omp_unset_lock(&vertexLocks[leafV.v]);
    }
    
    // ... tree 更新需要全局锁 ...
    
    for (auto &newLeaf : res.newLeaves) {
        #pragma omp critical(tree_update)
        {
            daf::Size newId = tree.addNode(newLeaf);
            // ...
        }
        
        for (auto i : tree.adj_list[newId]) {
            omp_set_lock(&vertexLocks[i.v]);
            treeGraphV.addNbr(i.v, {newId, i.isPivot});
            omp_unset_lock(&vertexLocks[i.v]);
        }
    }
}

// 清理锁
for (auto &lock : vertexLocks) omp_destroy_lock(&lock);
```

### 方案 3：优化 SDCT（最快达到 5×）

SDCT 占总时间的 78%，并行扩展性更好：
- 当前：2.4× 加速
- 潜力：10× 加速
- **如果 SDCT 达到 10×，总体可达 ~4.8× 加速**

## 调试工具

### 1. Valgrind（检测内存错误）
```bash
valgrind --leak-check=full --track-origins=yes \
  ./degeneracy_cliques graphs/com-dblp.edges 3 4
```

### 2. ThreadSanitizer（检测数据竞争）
```bash
# 编译时添加
-fsanitize=thread

# 运行
./degeneracy_cliques graphs/com-dblp.edges 3 4
```

### 3. GDB（调试崩溃）
```bash
gdb ./degeneracy_cliques
run graphs/com-dblp.edges 3 4
# 崩溃后
bt  # 查看调用栈
```

## 预期效果

### 如果批量操作成功
- tree 更新：串行（无法避免）
- treeGraphV 更新：并行（按顶点）
- **预期额外加速**：1.2-1.5×
- **总体加速**：1.69× × 1.3 ≈ **2.2×**

### 如果细粒度锁成功
- 完全并行更新
- **预期额外加速**：1.3-1.8×
- **总体加速**：1.69× × 1.5 ≈ **2.5×**

### 如果优化 SDCT
- SDCT：2.4× → 10×
- Nucleus：1.69×
- **总体加速**：**~4.8×**（接近 5× 目标！）

## 行动计划

### 立即执行（今天）
1. ✅ 实现方案 1 的步骤 1（简化版本）
2. ✅ 验证正确性
3. ✅ 逐步并行化
4. ✅ 测试性能

### 如果方案 1 失败
1. 尝试方案 2（细粒度锁）
2. 使用调试工具找出问题
3. 简化实现，确保稳定性

### 长期目标
1. 优化 SDCT（最有潜力）
2. 达到 5× 总体加速

## 总结

**我们已经找到了正确的方向！**

批量操作方案理论上是可行的，只是实现中遇到了问题。

通过：
1. 简化实现
2. 逐步测试
3. 使用调试工具
4. 仔细处理边界情况

**我们一定能够成功！**
