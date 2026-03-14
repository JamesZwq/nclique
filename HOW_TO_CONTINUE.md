# 如何继续优化实验

## 📊 当前状态

### 已完成
1. ✅ 创建了 16 个优化变体
2. ✅ 所有变体已同步到服务器
3. ✅ 实现了批量 treeGraphV 更新
4. ✅ 启用了 SDCT_Par
5. ✅ 创建了自动化测试框架

### 服务器上的文件
- `/home/wenqianz/nclique/optimization_variants/` - 16 个变体
- `/home/wenqianz/nclique/continuous_optimization.sh` - 持续优化脚本
- `/home/wenqianz/nclique/run_all_experiments_fixed.sh` - 完整测试脚本

## 🚀 如何手动运行实验

### 方法 1: 测试单个变体
```bash
ssh tods2
cd /home/wenqianz/nclique

# 备份当前文件
cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp \
   src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak

# 复制变体
cp optimization_variants/variant_01_schedule_dynamic_4.cpp \
   src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp

# 编译
cd build
make -j8 degeneracy_cliques

# 测试
cd /home/wenqianz/nclique
./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4

# 恢复
cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak \
   src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
```

### 方法 2: 运行自动化脚本
```bash
ssh tods2
cd /home/wenqianz/nclique

# 后台运行
nohup bash continuous_optimization.sh > continuous_opt.log 2>&1 &

# 监控进度
tail -f continuous_opt.log

# 查看结果
cat optimization_results/results.csv
```

### 方法 3: 批量测试所有变体
```bash
ssh tods2
cd /home/wenqianz/nclique

# 创建简单的循环脚本
cat > test_all.sh << 'SCRIPT'
#!/bin/bash
for variant in optimization_variants/variant_*.cpp; do
    name=$(basename $variant .cpp)
    echo "Testing $name..."
    
    cp $variant src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    cd build && make -j8 degeneracy_cliques > /dev/null 2>&1
    cd ..
    
    echo "$name," >> results.txt
    timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | \
        grep "Multi-threading" | tee -a results.txt
    
    echo "---" >> results.txt
done
SCRIPT

chmod +x test_all.sh
nohup bash test_all.sh &
```

## 📋 16 个优化变体说明

### Schedule 参数优化 (1-8)
1. **variant_01**: schedule(dynamic, 4) - 更小的 chunk
2. **variant_02**: schedule(dynamic, 32) - 更大的 chunk
3. **variant_03**: schedule(guided) - 自适应调度
4. **variant_04**: schedule(static) - 静态分配
5. **variant_05**: 批量 heap 更新
6. **variant_06**: 优化 Phase A schedule
7. **variant_07**: 同时优化两个 schedule
8. **variant_08**: Guided schedule for Phase B

### 高级优化 (9-16)
9. **variant_09**: Relaxed memory order
10. **variant_10**: 更大的 chunk (64)
11. **variant_11**: 更小的 chunk (8)
12. **variant_12**: Static with chunk (32)
13. **variant_13**: Phase A larger chunk (4)
14. **variant_14**: 两个阶段都用 guided
15. **variant_15**: 混合 schedule
16. **variant_16**: 优化 treeGraphV chunk (512)

## 🔍 如何分析结果

### 查看性能数据
```bash
ssh tods2
cd /home/wenqianz/nclique

# 查看所有结果
cat optimization_results/results.csv

# 找出最佳性能
grep "SUCCESS" optimization_results/results.csv | \
    sort -t',' -k3 -n | head -5

# 对比 baseline
grep "baseline" optimization_results/results.csv
```

### 性能改进计算
```bash
# 假设 baseline 是 436 ms
# 如果某个变体是 400 ms
# 改进 = (436 - 400) / 436 = 8.3%
```

## 🎯 预期结果

### 目标
找到能够带来 5-10% 改进的优化组合

### 如果找到好的变体
1. 记录变体名称和性能
2. 组合多个优化
3. 测试组合效果
4. 应用到最终代码

### 最终目标
- Nucleus: 1.69× → 2.0× (通过 peeling 优化)
- SDCT: 启用 SDCT_Par (8-10×)
- **总体: 5× 加速**

## 💡 下一步创建更多变体

### 如果需要更多优化
```bash
# 在本地
cd /Users/zhangwenqian/UNSW/pivoter

# 创建新变体
cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp \
   optimization_variants/variant_17_new_optimization.cpp

# 编辑文件实现优化
# ...

# 同步到服务器
scp optimization_variants/variant_17_new_optimization.cpp \
    tods2:/home/wenqianz/nclique/optimization_variants/

# 测试
ssh tods2 "cd /home/wenqianz/nclique && bash quick_test.sh variant_17_new_optimization"
```

## 📊 实时监控

### 创建监控脚本
```bash
# 在本地
cat > monitor.sh << 'SCRIPT'
#!/bin/bash
while true; do
    clear
    echo "=== Optimization Progress ==="
    echo "Time: $(date)"
    echo ""
    ssh tods2 "tail -20 /home/wenqianz/nclique/continuous_opt.log 2>/dev/null || \
               tail -20 /home/wenqianz/nclique/results.txt 2>/dev/null || \
               echo 'No logs yet'"
    echo ""
    echo "=== Results Summary ==="
    ssh tods2 "cat /home/wenqianz/nclique/optimization_results/results.csv 2>/dev/null | \
               tail -10 || echo 'No results yet'"
    sleep 30
done
SCRIPT
chmod +x monitor.sh
./monitor.sh
```

## 🔧 故障排除

### 如果实验失败
1. 检查 nCr.txt 文件是否存在
2. 确保从 /home/wenqianz/nclique 目录运行
3. 检查编译错误
4. 增加 timeout 时间

### 如果连接超时
1. 使用 nohup 后台运行
2. 定期检查日志文件
3. 使用更短的测试脚本

## 📝 记录模板

### 性能记录
```
Variant: variant_XX
Description: [优化说明]
com-dblp: XXX ms (baseline: 436 ms)
Improvement: XX%
Status: [保留/回滚]
Notes: [观察]
```

---

**所有工具和脚本已准备就绪，可以开始持续优化实验！**
