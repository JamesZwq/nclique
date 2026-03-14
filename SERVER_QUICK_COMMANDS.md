## 🚀 服务器测试 - 快速命令

### 方法1: 一键部署（推荐）

```bash
# 1. 上传代码到服务器
scp -r /Users/zhangwenqian/UNSW/pivoter tods2:~/

# 2. SSH到服务器并运行
ssh tods2 << 'EOF'
cd ~/pivoter
chmod +x deploy_server.sh test_server_32threads.sh
./deploy_server.sh
EOF
```

### 方法2: 手动步骤

```bash
# 1. 上传代码
scp -r /Users/zhangwenqian/UNSW/pivoter tods2:~/

# 2. SSH到服务器
ssh tods2

# 3. 进入目录
cd ~/pivoter

# 4. 设置权限
chmod +x deploy_server.sh test_server_32threads.sh

# 5. 运行部署脚本
./deploy_server.sh
```

### 方法3: 分步执行

```bash
# SSH到服务器
ssh tods2
cd ~/pivoter

# 编译
mkdir -p build-ultra && cd build-ultra
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc) test_ultra_parallel
cd ..

# 生成测试图
python3 generate_large_test_graphs.py

# 运行测试
chmod +x test_server_32threads.sh
./test_server_32threads.sh
```

### 快速测试单个配置

```bash
# 测试8线程
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 8

# 测试16线程
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 16

# 测试32线程
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 32
```

### 查看结果

```bash
# 查看最新结果
cat server_performance_results_*.txt | tail -50

# 查看加速比
grep "Speedup" server_performance_results_*.txt
```

### 预期输出

```
========================================
Testing: Medium (500 nodes)
========================================

Testing with 1 thread(s)...
  Reference (1 thread): 128ms
  Ultra Parallel (1 thread): 104ms (baseline)

Testing with 8 thread(s)...
  Ultra Parallel (8 threads): 39ms
    Speedup vs 1 thread: 2.66x
    Speedup vs Reference: 3.28x

Testing with 16 thread(s)...
  Ultra Parallel (16 threads): 25ms
    Speedup vs 1 thread: 4.16x
    Speedup vs Reference: 5.12x

Testing with 32 thread(s)...
  Ultra Parallel (32 threads): 18ms
    Speedup vs 1 thread: 5.77x
    Speedup vs Reference: 7.11x
```

### 故障排除

如果遇到问题：

```bash
# 检查编译日志
cat build-ultra/cmake_output.log
cat build-ultra/make_output.log

# 检查CPU核心数
nproc

# 检查OpenMP
echo $OMP_NUM_THREADS
export OMP_NUM_THREADS=32

# 重新编译
rm -rf build-ultra
./deploy_server.sh
```

### 成功标准

- ✓ 编译成功
- ✓ 正确性验证通过
- ✓ 8线程: >= 3.0x 加速
- ✓ 16线程: >= 4.0x 加速
- ✓ 32线程: >= 5.0x 加速

---

**准备就绪！** 所有脚本已配置完成，可以直接在服务器上运行。
