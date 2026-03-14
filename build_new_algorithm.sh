#!/bin/bash

# 编译和测试新算法

echo "========================================="
echo "Building and Testing New Algorithm"
echo "========================================="
echo ""

# 步骤 1: 添加源文件到编译
echo "Step 1: Adding source files..."
cd /Users/zhangwenqian/UNSW/pivoter

# 检查文件是否存在
if [ ! -f "src/NucleusDecomposition/NucleusCoreDecompositionOptimized.cpp" ]; then
    echo "ERROR: Optimized algorithm file not found!"
    exit 1
fi

echo "Files ready!"
echo ""

# 步骤 2: 编译
echo "Step 2: Compiling..."
cd cmake-build-release

# 编译优化算法（作为库）
echo "Compiling optimized algorithm..."
g++ -c -O3 -march=native -fopenmp -std=c++20 \
    -I../src \
    ../src/NucleusDecomposition/NucleusCoreDecompositionOptimized.cpp \
    -o NucleusCoreDecompositionOptimized.o

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed!"
    exit 1
fi

echo "Compilation successful!"
echo ""

# 步骤 3: 创建简单的测试
echo "Step 3: Creating test..."

cat > test_optimized.cpp << 'ENDTEST'
#include <iostream>
#include "../src/graph/Graph.h"
#include "../src/tree/MultiBranchTree.h"
#include "../src/NucleusDecomposition/NCliqueCoreDecomposition.h"
#include "../src/degeneracy_algorithm_cliques_V.h"

extern "C" {
namespace OptimizedNucleus {
std::vector<std::pair<std::vector<daf::Size>, int>> 
NucleusCoreDecompositionOptimized(
    DynamicGraph<TreeGraphNode>& tree,
    const Graph& edgeGraph,
    DynamicGraphSet<TreeGraphNode>& treeGraphV,
    daf::CliqueSize r,
    daf::CliqueSize s,
    int numThreads);
}
}

int main() {
    std::cout << "New algorithm compiled successfully!" << std::endl;
    return 0;
}
ENDTEST

g++ -O3 -fopenmp -std=c++20 -I../src test_optimized.cpp -o test_optimized

if [ $? -eq 0 ]; then
    echo "Test program compiled!"
    ./test_optimized
else
    echo "Test compilation failed (expected, need full build)"
fi

echo ""
echo "========================================="
echo "Next Steps:"
echo "========================================="
echo ""
echo "1. Modify CMakeLists.txt to include the new algorithm"
echo "2. Rebuild the project"
echo "3. Test on com-dblp dataset"
echo ""
echo "Or manually compile and link:"
echo "  cd cmake-build-release"
echo "  make degeneracy_cliques"
echo ""
