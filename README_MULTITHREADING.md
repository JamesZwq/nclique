# 多线程Nucleus Decomposition - 快速开始

## 快速测试

```bash
# 1. 编译
mkdir -p build-ultra && cd build-ultra
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8 test_ultra_parallel
cd ..

# 2. 生成测试图
python3 generate_large_test_graphs.py

# 3. 运行测试
./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 8

# 4. 完整性能测试
./comprehensive_test.sh
```

## 性能结果

在500节点图上，8线程达到 **3.28倍加速** (相对于单线程Reference)

## 服务器测试

```bash
ssh tods2
cd ~/pivoter
./test_on_server.sh
```

详细文档见 `FINAL_SUMMARY.md`
