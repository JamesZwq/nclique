#!/bin/bash

# 在服务器上运行大图测试

echo "=========================================="
echo "Running LARGE GRAPH tests on server"
echo "=========================================="
echo ""

# 1. 生成大图
echo "[1/2] Generating large test graphs..."
python3 generate_server_test_graphs.py

echo ""
echo "[2/2] Running performance tests..."
chmod +x test_server_large_graphs.sh
./test_server_large_graphs.sh

echo ""
echo "=========================================="
echo "Complete! Check server_large_results_*.txt"
echo "=========================================="
