#!/bin/bash

# 编译和测试Ultra Parallel算法

set -e

echo "=========================================="
echo "Building Ultra Parallel Algorithm"
echo "=========================================="

# 创建build目录
mkdir -p build-ultra
cd build-ultra

# 配置CMake（Release模式，启用优化）
cmake -DCMAKE_BUILD_TYPE=Release ..

# 编译
echo ""
echo "Compiling..."
make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

echo ""
echo "=========================================="
echo "Build Complete!"
echo "=========================================="

# 返回根目录
cd ..

# 检查是否有测试图
if [ ! -f "new_small_garph.edges" ]; then
    echo "Warning: new_small_garph.edges not found"
    echo "Please provide a test graph file"
    exit 1
fi

echo ""
echo "=========================================="
echo "Running Tests"
echo "=========================================="

# 测试不同线程数
for threads in 1 2 4 8; do
    echo ""
    echo "=========================================="
    echo "Testing with $threads threads"
    echo "=========================================="
    ./build-ultra/bin/test_ultra_parallel new_small_garph.edges 3 4 $threads
done

echo ""
echo "=========================================="
echo "All Tests Complete!"
echo "=========================================="
