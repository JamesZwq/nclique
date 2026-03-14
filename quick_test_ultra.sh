#!/bin/bash

# 简单的测试脚本

set -e

echo "=========================================="
echo "Testing Ultra Parallel Algorithm"
echo "=========================================="

# 使用小图测试
GRAPH="new_small_garph.edges"
R=3
S=4

echo ""
echo "Test 1: Single thread baseline"
echo "=========================================="
./build-ultra/bin/test_ultra_parallel $GRAPH $R $S 1

echo ""
echo ""
echo "Test 2: 2 threads"
echo "=========================================="
./build-ultra/bin/test_ultra_parallel $GRAPH $R $S 2

echo ""
echo ""
echo "Test 3: 4 threads"
echo "=========================================="
./build-ultra/bin/test_ultra_parallel $GRAPH $R $S 4

echo ""
echo ""
echo "Test 4: 8 threads"
echo "=========================================="
./build-ultra/bin/test_ultra_parallel $GRAPH $R $S 8

echo ""
echo "=========================================="
echo "All tests complete!"
echo "=========================================="
