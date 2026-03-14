#!/bin/bash
# 在服务器上编译和测试

cd ~/pivoter_new

echo "Compiling Advanced Parallel Algorithm..."
g++ -O3 -march=native -fopenmp -std=c++20 -I. -Isrc \
    test_advanced_parallel.cpp \
    src/NucleusDecomposition/NucleusCoreDecompositionAdvancedParallel.cpp \
    src/NucleusDecomposition/NucleusCoreDecompositionRemoveSclique.cpp \
    src/NucleusDecomposition/NCliqueCoreDecomposition.cpp \
    src/graph/graph.cpp \
    src/Global/Global.cpp \
    src/degeneracy_algorithm_cliques_V.cpp \
    src/degeneracy_helper.cpp \
    src/degeneracy_cliques.cpp \
    src/LinkedList.cpp \
    src/MemoryManager.cpp \
    src/misc.cpp \
    -o test_advanced 2>&1 | head -100

if [ $? -eq 0 ]; then
    echo "✓ Compilation successful!"
    echo ""
    echo "Running quick test with toy graph (4 threads)..."
    ./test_advanced toy.edges 3 4 4
else
    echo "✗ Compilation failed"
    exit 1
fi
