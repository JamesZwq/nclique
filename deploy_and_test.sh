#!/bin/bash

# 一键部署和测试脚本 - 用于服务器

echo "=========================================="
echo "Nucleus Decomposition Multi-threading"
echo "Deployment and Testing Script"
echo "=========================================="
echo ""

# 检查是否在正确的目录
if [ ! -f "CMakeLists.txt" ]; then
    echo "Error: CMakeLists.txt not found"
    echo "Please run this script from the project root directory"
    exit 1
fi

# 1. 清理旧的构建
echo "[1/6] Cleaning old build..."
rm -rf build-ultra
mkdir -p build-ultra

# 2. 配置CMake
echo "[2/6] Configuring CMake..."
cd build-ultra
cmake -DCMAKE_BUILD_TYPE=Release .. > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Error: CMake configuration failed"
    exit 1
fi

# 3. 编译
echo "[3/6] Compiling (this may take a while)..."
make -j8 test_ultra_parallel > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed"
    exit 1
fi
cd ..

echo "✓ Compilation successful"
echo ""

# 4. 生成测试图
echo "[4/6] Generating test graphs..."
if [ ! -f "test_medium_500.edges" ]; then
    python3 generate_large_test_graphs.py > /dev/null 2>&1
    echo "✓ Test graphs generated"
else
    echo "✓ Test graphs already exist"
fi
echo ""

# 5. 快速验证测试
echo "[5/6] Running quick verification test..."
OUTPUT=$(./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 1 2>&1)
if echo "$OUTPUT" | grep -q "Result sizes match"; then
    echo "✓ Correctness verification passed"
else
    echo "✗ Correctness verification failed"
    exit 1
fi
echo ""

# 6. 性能测试
echo "[6/6] Running performance tests..."
echo "=========================================="
echo ""

for T in 1 2 4 8; do
    echo "Testing with $T thread(s)..."
    ./build-ultra/bin/test_ultra_parallel test_medium_500.edges 3 4 $T 2>&1 | \
        grep -E "(Reference|Ultra Parallel|Speedup)" | head -3
    echo ""
done

echo "=========================================="
echo "Deployment and testing complete!"
echo ""
echo "Next steps:"
echo "  - Review FINAL_SUMMARY.md for detailed results"
echo "  - Run ./comprehensive_test.sh for detailed benchmarks"
echo "  - Test on larger graphs: test_large_1000.edges"
echo "=========================================="
