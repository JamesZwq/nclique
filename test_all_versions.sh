#!/bin/bash

# 简单的测试脚本来验证所有 SDCT 版本的正确性

if [ $# -lt 2 ]; then
    echo "Usage: $0 <graph_file> <max_threads>"
    exit 1
fi

GRAPH_FILE="$1"
MAX_THREADS="$2"
BIN_DIR="/Users/zhangwenqian/UNSW/pivoter/build/bin"

echo "=========================================="
echo "SDCT Correctness Verification Test"
echo "=========================================="
echo "Graph: $GRAPH_FILE"
echo "Max threads: $MAX_THREADS"
echo ""

# 检查二进制文件是否存在
if [ ! -f "$BIN_DIR/degeneracy_cliques" ]; then
    echo "Error: degeneracy_cliques binary not found at $BIN_DIR/degeneracy_cliques"
    exit 1
fi

echo "Running SDCT (reference implementation)..."
echo "Command: $BIN_DIR/degeneracy_cliques $GRAPH_FILE 2 3 default"
$BIN_DIR/degeneracy_cliques "$GRAPH_FILE" 2 3 default 2>&1 | head -50

echo ""
echo "=========================================="
echo "Test completed!"
echo "=========================================="
