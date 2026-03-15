#!/bin/bash

GRAPH="graphs/com-dblp.edges"
BIN="./build/bin/degeneracy_cliques"

echo "=========================================="
echo "Local Performance Test"
echo "=========================================="

# 测试 1 线程
echo "1 thread:"
export OMP_NUM_THREADS=1
time ./build/bin/degeneracy_cliques "$GRAPH" 2 3 2>&1 | grep "Tree Build"

echo ""
echo "4 threads:"
export OMP_NUM_THREADS=4
time ./build/bin/degeneracy_cliques "$GRAPH" 2 3 2>&1 | grep "Tree Build"

echo ""
echo "8 threads:"
export OMP_NUM_THREADS=8
time ./build/bin/degeneracy_cliques "$GRAPH" 2 3 2>&1 | grep "Tree Build"
