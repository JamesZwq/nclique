#!/bin/bash

# 性能测试脚本 - 测试不同线程数下各版本的性能

GRAPH="/data/wenqianz/com-dblp.edges"
BIN_DIR="/home/wenqianz/nclique/build/bin"
BINARY="$BIN_DIR/verify_all_sdct"

if [ ! -f "$BINARY" ]; then
    echo "Error: Binary not found at $BINARY"
    exit 1
fi

if [ ! -f "$GRAPH" ]; then
    echo "Error: Graph not found at $GRAPH"
    exit 1
fi

echo "=========================================="
echo "SDCT Performance Test"
echo "=========================================="
echo "Graph: $GRAPH"
echo "Binary: $BINARY"
echo ""

# 运行一次完整测试（用于验证正确性）
echo "Running full verification test..."
timeout 1800 "$BINARY" "$GRAPH" 2>&1 | tail -30

echo ""
echo "Test completed!"
