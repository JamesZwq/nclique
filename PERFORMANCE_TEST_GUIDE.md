# SDCT 性能测试脚本使用指南

## 脚本说明

本项目提供两个性能测试脚本：

### 1. `perf_test_server.sh` - 单版本性能测试

测试 SDCT_Par5 在不同线程数下的性能。

**使用方法：**

```bash
# 在服务器上运行
ssh tods2
cd /home/wenqianz/nclique_tmp

# 复制脚本（如果还没有）
cp /path/to/perf_test_server.sh .

# 运行脚本
chmod +x perf_test_server.sh
./perf_test_server.sh
```

**测试线程数：** 1, 2, 4, 8, 16, 32, 64

**输出示例：**
```
Threads | Time (ms) | Speedup
--------|-----------|--------
1       | 12345     | 1.00x
2       | 6500      | 1.90x
4       | 3400      | 3.63x
8       | 1800      | 6.86x
16      | 950       | 12.99x
32      | 520       | 23.74x
64      | 300       | 41.15x
```

### 2. `compare_versions.sh` - 版本对比测试

验证所有 SDCT 版本（Par2, Par3, Par4, Par5）的正确性。

**使用方法：**

```bash
# 在服务器上运行
ssh tods2
cd /home/wenqianz/nclique_tmp

# 复制脚本
cp /path/to/compare_versions.sh .

# 运行脚本
chmod +x compare_versions.sh
./compare_versions.sh
```

**输出示例：**
```
Testing SDCT...
SDCT clique counts: k1=317080 k2=1049866 ...

Testing SDCT_Par2...
SDCT_Par2 12 threads
Result: ✓ CORRECT

Testing SDCT_Par3...
Result: ✓ CORRECT

Testing SDCT_Par4...
Result: ✓ CORRECT

Testing SDCT_Par5...
Result: ✓ CORRECT

✓ ALL TESTS PASSED!
```

## 快速开始

### 步骤 1: 连接到服务器

```bash
ssh tods2
```

### 步骤 2: 进入项目目录

```bash
cd /home/wenqianz/nclique_tmp
```

### 步骤 3: 拉取最新代码

```bash
git pull origin main
```

### 步骤 4: 运行性能测试

```bash
# 方法 1: 测试单版本性能
bash perf_test_server.sh

# 方法 2: 验证所有版本正确性
bash compare_versions.sh
```

## 手动测试命令

如果脚本有问题，可以手动运行以下命令：

### 编译

```bash
cd /home/wenqianz/nclique_tmp
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make degeneracy_cliques -j16
make verify_all_sdct -j16
```

### 运行验证

```bash
./bin/verify_all_sdct /data/wenqianz/com-dblp.edges
```

### 运行性能测试（不同线程数）

```bash
# 单线程
export OMP_NUM_THREADS=1
time ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 2 3

# 8 线程
export OMP_NUM_THREADS=8
time ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 2 3

# 16 线程
export OMP_NUM_THREADS=16
time ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 2 3

# 32 线程
export OMP_NUM_THREADS=32
time ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 2 3

# 64 线程
export OMP_NUM_THREADS=64
time ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 2 3
```

## 预期结果

### 正确性
- ✓ 所有版本（Par2, Par3, Par4, Par5）都应该通过验证
- ✓ cliqueCount 应该完全一致（相对误差 < 1e-10）

### 性能
- 单线程基准：~12-15 秒
- 64 线程加速比：40-50x
- 线性加速比应该随线程数增加而增加

## 故障排除

### 问题 1: CMake 配置超时

**解决方案：** 使用已编译的二进制文件，或在后台运行编译

```bash
nohup bash -c 'cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j16' > /tmp/build.log 2>&1 &
```

### 问题 2: 编译失败

**解决方案：** 清理并重新编译

```bash
cd /home/wenqianz/nclique_tmp
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j16
```

### 问题 3: 脚本权限不足

**解决方案：** 添加执行权限

```bash
chmod +x perf_test_server.sh
chmod +x compare_versions.sh
```

## 输出文件

- `/tmp/verify_results.txt` - 版本验证结果
- `/tmp/sdct_perf_results.txt` - 性能测试结果

## 注意事项

1. **图文件位置：** `/data/wenqianz/com-dblp.edges`
2. **项目位置：** `/home/wenqianz/nclique_tmp`
3. **线程数：** 服务器有 96 个核心，可以测试到 64 线程
4. **运行时间：** 完整测试可能需要 30-60 分钟

## 联系方式

如有问题，请检查：
1. 图文件是否存在
2. 项目代码是否最新（git pull）
3. 编译是否成功
4. 线程数设置是否正确
