# 持续优化进度报告

## 目标：达到 5× 加速

当前最佳结果：web-Google 32线程 **1.63×**（距离目标 33%）

---

## 已尝试的优化

### 1. 优化 countingPerRClique ✅
- 移除 `#pragma omp critical`
- 使用线程局部数组 + 并行 reduction
- **效果**：1.68-2.72× 加速

### 2. 优化 Part A（收集 leaves）✅
- 移除 `#pragma omp critical`
- 使用线程局部缓冲区
- **效果**：有一定改善

### 3. 并行更新 support（std::atomic）✅
- 使用 `std::atomic<double>` + `compare_exchange_weak`
- **效果**：
  - web-Google 32线程：1.07× → **1.63×**（提升 52%）
  - com-dblp：性能下降（atomic 开销太大）

### 4. 并行更新 support（omp atomic）❌
- 使用 `#pragma omp atomic`
- **效果**：比 std::atomic 更慢
  - web-Google 16线程：1.08×（比 1.63× 差很多）

### 5. 完全并行化 tree/treeGraphV 更新 ❌
- 使用多个 critical sections
- **效果**：程序崩溃（死锁或数据竞争）

---

## 当前瓶颈分析

### 时间分布（web-Google，单线程）
- countingPerRClique: ~67 ms（10%）✅ 已优化
- 主循环: ~650 ms（90%）
  - Part A（收集 leaves）: ~50 ms ✅ 已优化
  - Part B（处理 leaves）: ~200 ms ✅ 已并行
  - **Update tree/treeGraphV: ~300 ms** ❌ 串行瓶颈
  - **Update support: ~100 ms** ⚠️ 部分优化

### 最大的瓶颈：tree/treeGraphV 更新（~300 ms，46%）

这部分完全串行，包括：
1. `treeGraphV.removeNbr()` - 移除旧 leaf
2. `tree.addNode()` - 添加新 leaf
3. `treeGraphV.addNbr()` - 添加新 leaf 的连接
4. `tree.removeNode()` - 移除旧 leaf

---

## 下一步优化策略

### 策略 A：细粒度锁（推荐）
为 tree 和 treeGraphV 的每个节点/边创建细粒度锁

**优点**：
- 不同 leaves 的更新可以并行
- 只有访问相同节点时才会冲突

**缺点**：
- 实现复杂
- 锁的开销

**预期加速**：2-3×

### 策略 B：无锁数据结构
使用无锁的并发数据结构

**优点**：
- 无锁，性能更好
- 完全并行

**缺点**：
- 实现非常复杂
- 需要重写 tree 和 treeGraphV

**预期加速**：3-5×

### 策略 C：批量更新
收集所有更新，然后批量应用

**优点**：
- 减少同步次数
- 可以优化更新顺序

**缺点**：
- 需要额外内存
- 实现复杂

**预期加速**：1.5-2×

### 策略 D：分区并行
将 tree/treeGraphV 分成多个分区，每个分区独立更新

**优点**：
- 减少冲突
- 易于实现

**缺点**：
- 负载不均衡
- 跨分区的更新仍需同步

**预期加速**：1.5-2×

---

## 实现计划

### 第一步：实现策略 D（分区并行）
最简单，风险最小

### 第二步：如果不够，实现策略 A（细粒度锁）
更复杂，但预期效果更好

### 第三步：如果还不够，考虑策略 B（无锁数据结构）
最复杂，但可能达到 5× 目标

---

## 当前状态

- ✅ countingPerRClique: 已优化
- ✅ Part A: 已优化
- ✅ Part B: 已并行
- ⚠️ Update support: 部分优化（大图有效）
- ❌ Update tree/treeGraphV: **最大瓶颈，完全串行**

**下一个目标**：并行化 tree/treeGraphV 更新，目标加速 2-3×
