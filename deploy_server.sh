#!/bin/bash

# 服务器部署和测试脚本 - 完整版

echo "=========================================="
echo "Server Deployment and Testing"
echo "Nucleus Decomposition Multi-threading"
echo "=========================================="
echo ""

# 1. 环境检查
echo "[1/7] Checking environment..."
echo "  Hostname: $(hostname)"
echo "  CPU Cores: $(nproc)"
echo "  Available memory: $(free -h | grep Mem | awk '{print $2}')"
echo ""

# 检查必要的工具
if ! command -v cmake &> /dev/null; then
    echo "Error: cmake not found. Please install cmake."
    exit 1
fi

if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found. Please install python3."
    exit 1
fi

echo "✓ Environment check passed"
echo ""

# 2. 清理旧构建
echo "[2/7] Cleaning old build..."
rm -rf build-ultra
mkdir -p build-ultra

# 3. 配置CMake
echo "[3/7] Configuring CMake..."
cd build-ultra
cmake -DCMAKE_BUILD_TYPE=Release .. > cmake_output.log 2>&1
if [ $? -ne 0 ]; then
    echo "Error: CMake configuration failed"
    echo "See cmake_output.log for details"
    exit 1
fi
cd ..
echo "✓ CMake configured"
echo ""

# 4. 编译
echo "[4/7] Compiling (this may take a few minutes)..."
cd build-ultra
make -j$(nproc) test_ultra_parallel > make_output.log 2>&1
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed"
    echo "See make_output.log for details"
    exit 1
fi
cd ..
echo "✓ Compilation successful"
echo ""

# 5. 生成测试图
echo "[5/7] Generating test graphs..."
if [ ! -f "test_medium_500.edges" ]; then
    python3 generate_large_test_graphs.py > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate test graphs"
        exit 1
    fi
    echo "✓ Test graphs generated"
else
    echo "✓ Test graphs already exist"
fi
echo ""

# 6. 快速验证
echo "[6/7] Running correctness verification..."
OUTPUT=$(./build-ultra/bin/test_ultra_parallel test_small.edges 3 4 1 2>&1)
if echo "$OUTPUT" | grep -q "Result sizes match"; then
    echo "✓ Correctness verification passed"
else
    echo "✗ Correctness verification failed"
    exit 1
fi
echo ""

# 7. 性能测试
echo "[7/7] Running performance tests..."
echo "=========================================="
echo ""

./test_server_32threads.sh

echo ""
echo "=========================================="
echo "Deployment and testing complete!"
echo ""
echo "Next steps:"
echo "  - Review the results file: server_performance_results_*.txt"
echo "  - Check FINAL_CONFIRMATION.md for detailed analysis"
echo "  - Run additional tests if needed"
echo "=========================================="
