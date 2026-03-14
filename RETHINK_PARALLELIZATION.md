# 重新思考：Nucleus Decomposition 的并行化

## 关键洞察：我可能理解错了

让我重新审视论文中的算法描述...

### 当前理解的问题

我一直认为 tree/treeGraphV 的更新必须串行，因为：
1. 不同 leaves 可能访问相同的节点
2. 会有数据竞争

**但是，等等！让我重新思考...**

### 新的洞察

#### 观察 1：不同 leaf 的独立性

论文中提到：
- 每个 leaf 对应一个 (r,s)-clique tree 的节点
- 移除一个 r-clique 后，只影响**包含它的 leaves**
- 不同的 leaves 在 tree 中是**独立的节点**

**关键问题**：不同 leaves 的更新真的会冲突吗？

#### 观察 2：treeGraphV 的结构

`treeGraphV` 是一个从 vertex 到 leaf 的映射：
- `treeGraphV[v]` = 包含顶点 v 的所有 leaves

当我们：
1. 移除 leaf L1，从 treeGraphV[v] 中删除 L1
2. 添加 leaf L2，向 treeGraphV[v] 中添加 L2

**关键问题**：如果两个线程同时操作不同的 leaves，它们会访问相同的 treeGraphV[v] 吗？

#### 观察 3：可能的并行策略

**策略 A：按顶点分区**
- 将顶点分成多个分区
- 每个线程负责一个分区的 treeGraphV 更新
- 只有当 leaf 包含该分区的顶点时，才由该线程处理

**策略 B：延迟合并**
- 每个线程维护自己的 tree 和 treeGraphV 副本
- 并行处理所有 leaves
- 最后合并所有副本

**策略 C：读写分离**
- 读取阶段：并行读取所有需要的数据
- 写入阶段：使用细粒度锁或原子操作

### 让我重新检查论文

论文中是否提到了并行化的方法？
让我仔细阅读...
