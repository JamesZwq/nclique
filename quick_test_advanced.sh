#!/bin/bash

# 快速测试脚本 - 在本地或服务器上快速验证编译和基本功能

set -e

echo "Quick test of Advanced Parallel Algorithm"
echo ""

# 使用最小的测试数据
TEST_GRAPH="toy.edges"
R=3
S=4
THREADS=4

if [ ! -f "$TEST_GRAPH" ]; then
    echo "ERROR: Test graph $TEST_GRAPH not found"
    exit 1
fi

echo "Compiling..."
g++ -O3 -march=native -fopenmp -std=c++20 \
    -I. \
    test_advanced_parallel.cpp \
    src/NucleusDecomposition/NucleusCoreDecompositionAdvancedParallel.cpp \
    src/NucleusDecomposition/NucleusCoreDecompositionRemoveSclique.cpp \
    src/NucleusDecomposition/NCliqueCoreDecomposition.cpp \
    src/graph/Graph.cpp \
    src/Global/Global.cpp \
    src/degeneracy_algorithm_cliques_V.cpp \
    src/degeneracy_helper.cpp \
    src/degeneracy_cliques.cpp \
    src/LinkedList.cpp \
    src/MemoryManager.cpp \
    src/misc.cpp \
    -o test_advanced_quick 2>&1 | head -20

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful"
    echo ""
    echo "Running quick test..."
    ./test_advanced_quick "$TEST_GRAPH" "$R" "$S" "$THREADS"
else
    echo "✗ Compilation failed"
    exit 1
fi
