## 快速开始指南

### 本地测试（已完成）

✅ 已在本地完成测试，达到 **3.28倍加速** (8线程 vs Reference单线程)

### 服务器测试步骤

#### 方法1: 一键部署（推荐）

```bash
# 1. 上传整个项目到服务器
scp -r /Users/zhangwenqian/UNSW/pivoter tods2:~/

# 2. SSH到服务器
ssh tods2

# 3. 运行一键部署脚本
cd ~/pivoter
./deploy_and_test.sh
```

#### 方法2: 手动步骤

```bash
# 1. SSH到服务器
ssh tods2

# 2. 进入项目目录
cd ~/pivoter

# 3. 编译
mkdir -p build-ultra && cd build-ultra
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8 test_ultra_parallel
cd ..

# 4. 生成测试图
python3 generate_large_test_graphs.py

# 5. 运行测试
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 8

# 6. 完整性能测试
./comprehensive_test.sh
```

### 预期结果

在服务器上应该能看到：
- ✓ 编译成功
- ✓ 正确性验证通过
- ✓ 8线程加速比 >= 3.0x

### 关键文件

- `FINAL_SUMMARY.md` - 完整工作总结
- `PERFORMANCE_REPORT.md` - 性能测试报告
- `deploy_and_test.sh` - 一键部署脚本
- `comprehensive_test.sh` - 详细性能测试

### 故障排除

如果遇到问题：

1. **编译错误**: 检查Boost和OpenMP是否安装
2. **运行错误**: 确保测试图已生成
3. **性能不佳**: 检查CPU核心数和OpenMP设置

```bash
# 检查OpenMP
echo $OMP_NUM_THREADS

# 设置线程数
export OMP_NUM_THREADS=8
```

### 联系

如有问题，查看 `FINAL_SUMMARY.md` 中的详细文档。
