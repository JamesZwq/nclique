## 服务器测试指南

### 快速开始

```bash
# 1. 上传代码到服务器
scp -r /Users/zhangwenqian/UNSW/pivoter tods2:~/

# 2. SSH到服务器
ssh tods2

# 3. 进入项目目录
cd ~/pivoter

# 4. 运行一键部署脚本
chmod +x deploy_server.sh test_server_32threads.sh
./deploy_server.sh
```

### 测试配置

- **线程数**: 1, 2, 4, 8, 16, 32
- **测试图**: 
  - Small (50 nodes)
  - Medium (500 nodes)
  - Large (1000 nodes)
- **参数**: r=3, s=4
- **运行次数**: 每个配置运行3次取平均

### 预期结果

在服务器上（更好的CPU和内存带宽），预期能看到：

- ✓ 8线程: 3-5倍加速
- ✓ 16线程: 4-7倍加速
- ✓ 32线程: 5-10倍加速（如果有足够的核心）

### 输出文件

测试完成后会生成：
- `server_performance_results_YYYYMMDD_HHMMSS.txt` - 详细结果

### 手动测试

如果需要手动测试特定配置：

```bash
# 测试特定线程数
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 16

# 详细验证
./build-ultra/bin/test_detailed_verification test_medium_500.edges 3 4 32
```

### 故障排除

#### 编译错误

```bash
# 检查Boost
ls ~/bin/boost/include

# 检查编译器
g++ --version

# 查看详细错误
cat build-ultra/cmake_output.log
cat build-ultra/make_output.log
```

#### 运行错误

```bash
# 检查OpenMP
echo $OMP_NUM_THREADS

# 设置线程数
export OMP_NUM_THREADS=32

# 检查测试图
ls -lh test_*.edges
```

#### 性能不佳

```bash
# 检查CPU信息
lscpu

# 检查负载
top

# 检查内存
free -h
```

### 性能分析

如果需要详细的性能分析：

```bash
# 使用perf（如果可用）
perf stat -e cycles,instructions,cache-misses \
  ./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 32

# 使用time
time ./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 32
```

### 预期时间

- 编译: 2-5分钟
- 测试图生成: 10-30秒
- 完整测试: 5-15分钟（取决于图大小和线程数）

### 结果解读

查看生成的结果文件：

```bash
cat server_performance_results_*.txt
```

关注以下指标：
- **Speedup vs Reference**: 相对于原始单线程的加速比
- **Speedup vs 1 thread**: 相对于优化后单线程的加速比
- **效率**: 加速比 / 线程数（理想值接近1.0）

### 成功标准

- ✓ 编译成功，无错误
- ✓ 正确性验证通过
- ✓ 8线程加速比 >= 3.0x
- ✓ 16线程加速比 >= 4.0x
- ✓ 32线程加速比 >= 5.0x（如果有足够核心）

### 联系

如有问题，查看以下文档：
- `FINAL_CONFIRMATION.md` - 最终确认报告
- `CORRECTNESS_VERIFICATION_REPORT.md` - 正确性验证
- `PROJECT_COMPLETION_REPORT.md` - 完整技术文档
