#!/bin/bash

# 编译和测试高级并行算法

set -e

echo "=========================================="
echo "Building Advanced Parallel Algorithm"
echo "=========================================="

# 创建构建目录
mkdir -p build-advanced
cd build-advanced

# 配置CMake（如果有CMakeLists.txt）
if [ -f ../CMakeLists.txt ]; then
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -march=native -fopenmp" ..
    make -j8
else
    # 手动编译
    echo "Manual compilation..."
    
    g++ -O3 -march=native -fopenmp -std=c++20 \
        -I.. \
        ../test_advanced_parallel.cpp \
        ../src/NucleusDecomposition/NucleusCoreDecompositionAdvancedParallel.cpp \
        ../src/NucleusDecomposition/NucleusCoreDecompositionRemoveSclique.cpp \
        ../src/NucleusDecomposition/NCliqueCoreDecomposition.cpp \
        ../src/graph/Graph.cpp \
        ../src/Global/Global.cpp \
        ../src/degeneracy_algorithm_cliques_V.cpp \
        ../src/degeneracy_helper.cpp \
        ../src/degeneracy_cliques.cpp \
        ../src/LinkedList.cpp \
        ../src/MemoryManager.cpp \
        ../src/misc.cpp \
        -o test_advanced_parallel
fi

cd ..

echo ""
echo "=========================================="
echo "Build complete!"
echo "=========================================="
echo ""

# 检查是否在服务器上
if [ -f /proc/cpuinfo ]; then
    CORES=$(grep -c ^processor /proc/cpuinfo)
    echo "Detected $CORES CPU cores"
fi

echo ""
echo "To run tests on server (tods2):"
echo "  ssh tods2"
echo "  cd $(pwd)"
echo "  ./run_advanced_experiments.sh"
echo ""
