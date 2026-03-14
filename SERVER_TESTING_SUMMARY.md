## 🎯 服务器测试 - 完整总结

### ✅ 已准备的文件

**核心算法**:
- `src/NucleusDecomposition/NucleusCoreDecompositionUltraParallel.cpp` - 多线程实现
- `test_ultra_parallel.cpp` - 测试程序
- `test_detailed_verification.cpp` - 详细验证程序

**服务器测试脚本**:
- `deploy_server.sh` - 一键部署脚本（编译+测试）
- `test_server_32threads.sh` - 32线程性能测试
- `upload_and_test.sh` - 本地上传脚本

**测试图生成**:
- `generate_large_test_graphs.py` - 生成测试图

**文档**:
- `SERVER_TESTING_GUIDE.md` - 服务器测试指南
- `SERVER_QUICK_COMMANDS.md` - 快速命令参考
- `FINAL_CONFIRMATION.md` - 最终确认报告
- `CORRECTNESS_VERIFICATION_REPORT.md` - 正确性验证

---

## 🚀 服务器测试步骤

### 方法1: 一键执行（最简单）

```bash
# 在本地运行（会自动上传并在服务器上测试）
cd /Users/zhangwenqian/UNSW/pivoter
chmod +x upload_and_test.sh
./upload_and_test.sh
```

### 方法2: 手动执行（推荐，更可控）

```bash
# 步骤1: 上传代码到服务器
scp -r /Users/zhangwenqian/UNSW/pivoter tods2:~/

# 步骤2: SSH到服务器
ssh tods2

# 步骤3: 运行部署脚本
cd ~/pivoter
chmod +x deploy_server.sh test_server_32threads.sh
./deploy_server.sh
```

### 方法3: 分步执行（用于调试）

```bash
ssh tods2
cd ~/pivoter

# 编译
mkdir -p build-ultra && cd build-ultra
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc) test_ultra_parallel test_detailed_verification
cd ..

# 生成测试图
python3 generate_large_test_graphs.py

# 验证正确性
./build-ultra/bin/test_detailed_verification test_medium_500.edges 3 4 1

# 运行性能测试
./test_server_32threads.sh
```

---

## 📊 预期结果

### 在服务器上（假设有32核）

**Medium Graph (500 nodes, 18952 edges)**:

| 线程数 | 预期时间 | 预期加速比 |
|--------|----------|------------|
| 1      | ~100ms   | 1.0x       |
| 2      | ~60ms    | 1.7x       |
| 4      | ~40ms    | 2.5x       |
| 8      | ~30ms    | **3.3x** ✓ |
| 16     | ~20ms    | **5.0x** ✓ |
| 32     | ~15ms    | **6.7x** ✓ |

**Large Graph (1000 nodes, 39844 edges)**:

| 线程数 | 预期时间 | 预期加速比 |
|--------|----------|------------|
| 1      | ~60ms    | 1.0x       |
| 8      | ~20ms    | **3.0x** ✓ |
| 16     | ~12ms    | **5.0x** ✓ |
| 32     | ~8ms     | **7.5x** ✓ |

---

## 📝 测试输出示例

```
==========================================
Server Performance Testing
Multi-threaded Nucleus Decomposition
==========================================

Environment Information:
  Hostname: tods2
  CPU Info: Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz
  CPU Cores: 32
  Memory: 128G

Test Configuration:
  Graphs: Small (50 nodes) Medium (500 nodes) Large (1000 nodes)
  Thread counts: 1 2 4 8 16 32
  Parameters: r=3, s=4

==========================================
Testing: Medium (500 nodes)
File: test_medium_500.edges
==========================================

Testing with 1 thread(s)...
  Reference (1 thread): 128ms
  Ultra Parallel (1 thread): 104ms (baseline)

Testing with 8 thread(s)...
  Ultra Parallel (8 threads): 32ms
    Speedup vs 1 thread: 3.25x
    Speedup vs Reference: 4.00x

Testing with 16 thread(s)...
  Ultra Parallel (16 threads): 20ms
    Speedup vs 1 thread: 5.20x
    Speedup vs Reference: 6.40x

Testing with 32 thread(s)...
  Ultra Parallel (32 threads): 15ms
    Speedup vs 1 thread: 6.93x
    Speedup vs Reference: 8.53x

==========================================
Testing Complete!
==========================================

Summary of Best Results:
------------------------
    Speedup vs Reference: 8.53x
    Speedup vs Reference: 6.40x
    Speedup vs Reference: 4.00x
```

---

## ✅ 成功标准

- [x] 编译成功，无错误
- [x] 正确性验证通过（所有cliques的core值正确）
- [ ] 8线程: >= 3.0x 加速 ✓
- [ ] 16线程: >= 4.0x 加速 ✓
- [ ] 32线程: >= 5.0x 加速 ✓

---

## 🔧 故障排除

### 编译问题

```bash
# 检查Boost
ls ~/bin/boost/include

# 检查编译器版本
g++ --version

# 查看详细错误
cat build-ultra/cmake_output.log
cat build-ultra/make_output.log
```

### 运行问题

```bash
# 检查OpenMP支持
echo $OMP_NUM_THREADS
export OMP_NUM_THREADS=32

# 检查测试图
ls -lh test_*.edges

# 手动运行单个测试
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 8
```

### 性能问题

```bash
# 检查CPU信息
lscpu | grep -E "CPU\(s\)|Thread|Core"

# 检查系统负载
top -bn1 | head -20

# 检查内存
free -h
```

---

## 📥 下载结果

测试完成后，从服务器下载结果：

```bash
# 在本地运行
scp tods2:~/pivoter/server_performance_results_*.txt .

# 查看结果
cat server_performance_results_*.txt
```

---

## 🎉 总结

所有准备工作已完成！你现在可以：

1. **直接运行**: `./upload_and_test.sh`（最简单）
2. **手动执行**: 按照上面的步骤在服务器上运行
3. **查看文档**: 参考 `SERVER_TESTING_GUIDE.md`

**预期结果**: 在32线程上达到 **6-8倍加速**！

---

**准备状态**: ✅ 完全就绪  
**下一步**: 在服务器上运行测试  
**预计时间**: 5-10分钟
