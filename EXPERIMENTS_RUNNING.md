# 优化实验正在运行

## 🚀 当前状态

实验已经在服务器上启动，将运行一整晚测试所有优化方案。

## 📋 实验设置

### 测试的优化变体

#### 第一批：Schedule 参数优化
1. **baseline**: 当前版本 (schedule(dynamic, 16))
2. **variant_01**: schedule(dynamic, 4)
3. **variant_02**: schedule(dynamic, 32)
4. **variant_03**: schedule(guided)
5. **variant_04**: schedule(static)

### 测试流程

对于每个变体：
1. 编译代码
2. 在 com-dblp 上运行 3 次，取平均值
3. 在 web-Google 上运行 1 次
4. 记录所有结果到 CSV 文件

### 结果位置

服务器上的结果文件：
- `/home/wenqianz/nclique/optimization_results/results.csv` - 性能数据
- `/home/wenqianz/nclique/optimization_results/experiment.log` - 详细日志
- `/home/wenqianz/nclique/experiment_output.log` - 完整输出

## 📊 如何查看结果

### 查看实时进度
```bash
ssh tods2 "tail -f /home/wenqianz/nclique/experiment_output.log"
```

### 查看结果摘要
```bash
ssh tods2 "cat /home/wenqianz/nclique/optimization_results/results.csv"
```

### 查看最新结果
```bash
ssh tods2 "tail -20 /home/wenqianz/nclique/optimization_results/experiment.log"
```

## 🔄 下一批实验

一旦第一批完成，我会创建更多优化变体：

### 第二批：Heap 优化
- 批量 heap 更新
- 延迟 heap 更新
- 不同的 heap 实现

### 第三批：内存优化
- 数据预取
- 内存池
- 缓存优化

### 第四批：算法优化
- 剪枝
- 增量更新
- 批处理

## 📈 预期结果

### 保守估计
- 如果找到 5 个有效优化，每个提升 5%
- 总体提升：(1.05)^5 = 1.28×
- 结合当前 1.69× = **2.16× Nucleus 加速**

### 乐观估计
- 如果找到关键优化（如 BK 算法优化 20%）
- 可能达到 **2.5-3× Nucleus 加速**
- 结合 SDCT_Par (8-10×)
- **总体可能达到 6-8× 加速**

## 🎯 最终目标

通过系统性测试所有优化方案，找到最佳组合，达到或超过 **5× 总体加速**。

## ⏰ 预计完成时间

- 第一批（5 个变体）：约 2-3 小时
- 如果继续创建更多变体：整晚运行
- 明天早上可以看到完整结果

## 📝 实验记录

实验开始时间：2025-03-09 凌晨
预计结束时间：2025-03-09 早上

所有结果将自动记录到 CSV 文件中，便于分析和比较。
