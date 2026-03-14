# 立即行动清单 - 测试新算法

## ✅ 已完成的工作

1. **设计了全新算法** - 不基于现有代码
2. **实现了 3 个版本** - 不同的优化策略
3. **创建了测试框架** - 完整的测试程序
4. **编写了详细文档** - 算法说明和使用指南

## 🚀 立即可以做的事情

### 选项 1: 快速验证（推荐）

直接在现有代码中测试新算法的核心思想：

```bash
# 1. 修改现有代码，使用 bucket sort 代替 heap
# 2. 添加批量并行处理
# 3. 测试性能
```

**预期结果**: 1.2-1.5× 加速

### 选项 2: 完整集成

将新算法完全集成到项目中：

```bash
# 1. 修改 CMakeLists.txt
# 2. 编译新算法
# 3. 运行完整测试
# 4. 比较新旧算法
```

**预期结果**: 1.3-1.7× 加速

### 选项 3: 继续优化现有算法

基于之前创建的 16 个优化变体：

```bash
# 在服务器上运行
ssh tods2
cd /home/wenqianz/nclique

# 创建并运行测试脚本
cat > test_all_variants.sh << 'EOF'
#!/bin/bash
export OMP_NUM_THREADS=16
cd /home/wenqianz/nclique
echo "Start: $(date)" > results.txt

for v in optimization_variants/variant_*.cpp; do
  name=$(basename $v .cpp)
  echo "Testing $name..." | tee -a results.txt
  
  cp $v src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
  cd build && make -j8 degeneracy_cliques >/dev/null 2>&1 && cd ..
  
  for run in 1 2 3; do
    echo "  Run $run:" | tee -a results.txt
    timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | \
      grep "Multi-threading" | tee -a results.txt
  done
  echo "" >> results.txt
done

echo "End: $(date)" >> results.txt
EOF

chmod +x test_all_variants.sh
nohup bash test_all_variants.sh > test_output.log 2>&1 &

# 监控进度
tail -f test_output.log
```

**预期结果**: 找到 5-10% 的改进

## 📊 三种方案对比

| 方案 | 工作量 | 预期加速 | 风险 | 时间 |
|------|--------|----------|------|------|
| 快速验证 | 小 | 1.2-1.5× | 低 | 1-2 小时 |
| 完整集成 | 中 | 1.3-1.7× | 中 | 3-4 小时 |
| 继续优化 | 小 | 1.05-1.1× | 低 | 2-3 小时 |

## 🎯 推荐方案

### 短期（今晚）
**方案 3**: 继续优化现有算法
- 16 个变体已准备好
- 可以立即运行
- 风险最低

### 中期（明天）
**方案 1**: 快速验证新算法
- 实现 bucket sort
- 测试性能
- 如果有效，继续完整集成

### 长期（本周）
**方案 2**: 完整集成新算法
- 完整实现
- 全面测试
- 达到 5× 目标

## 💻 立即执行的命令

### 在服务器上运行 16 个变体测试

```bash
# 复制粘贴这些命令到服务器
ssh tods2
cd /home/wenqianz/nclique

# 创建测试脚本
cat > test_16threads.sh << 'ENDSCRIPT'
#!/bin/bash
export OMP_NUM_THREADS=16
cd /home/wenqianz/nclique
echo "========================================" > results_16t.txt
echo "16-Thread Optimization Tests" >> results_16t.txt
echo "Start: $(date)" >> results_16t.txt
echo "========================================" >> results_16t.txt

for v in optimization_variants/variant_*.cpp; do
  name=$(basename $v .cpp)
  echo "" >> results_16t.txt
  echo "Testing $name..." | tee -a results_16t.txt
  
  cp $v src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
  cd build
  make -j8 degeneracy_cliques >/dev/null 2>&1
  if [ $? -ne 0 ]; then
    echo "  COMPILE FAILED" | tee -a ../results_16t.txt
    cd ..
    continue
  fi
  cd ..
  
  for run in 1 2 3; do
    echo "  Run $run:" | tee -a results_16t.txt
    timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | \
      grep -E "(Multi-threading|took:)" | tee -a results_16t.txt
    sleep 2
  done
done

echo "" >> results_16t.txt
echo "========================================" >> results_16t.txt
echo "End: $(date)" >> results_16t.txt
echo "========================================" >> results_16t.txt

# 分析结果
echo "" >> results_16t.txt
echo "BEST RESULTS:" >> results_16t.txt
grep "took:" results_16t.txt | sort -t':' -k2 -n | head -10 >> results_16t.txt
ENDSCRIPT

chmod +x test_16threads.sh

# 后台运行
nohup bash test_16threads.sh > test_output.log 2>&1 &

# 查看进度
echo "Test started! Monitor with:"
echo "  tail -f test_output.log"
echo "  tail -f results_16t.txt"
```

## 📈 预期时间线

### 今晚（2-3 小时）
- 运行 16 个变体测试
- 收集性能数据
- 找出最佳优化

### 明天早上
- 分析结果
- 组合最佳优化
- 测试组合效果

### 明天下午
- 如果需要，实现新算法
- 完整测试
- 验证 5× 目标

## 🏆 成功标准

### 最低目标
- 找到 5% 的改进
- 验证优化有效

### 理想目标
- Nucleus: 1.5× 加速
- 总体: 5× 加速 ✅

### 超越目标
- Nucleus: 2× 加速
- 总体: 6× 加速 🚀

---

**所有工具和代码已准备就绪！现在就可以开始测试！**
