#!/bin/bash

# 服务器测试脚本 - 在tods2上运行

echo "=========================================="
echo "Server Performance Testing (tods2)"
echo "=========================================="
echo ""

# 检查是否在服务器上
if [ ! -f "/etc/hostname" ] || ! grep -q "tods2" /etc/hostname 2>/dev/null; then
    echo "Warning: This script is designed for tods2 server"
    echo "Current hostname: $(hostname)"
    echo ""
fi

# 设置线程数
export OMP_NUM_THREADS=8

# 测试配置
GRAPHS=("test_medium_500.edges" "test_large_1000.edges")
R=3
S=4

echo "Configuration:"
echo "  Max threads: 8"
echo "  Graphs: ${GRAPHS[@]}"
echo "  Parameters: r=$R, s=$S"
echo ""

for GRAPH in "${GRAPHS[@]}"; do
    if [ ! -f "$GRAPH" ]; then
        echo "Error: $GRAPH not found"
        continue
    fi
    
    echo "=========================================="
    echo "Testing: $GRAPH"
    echo "=========================================="
    echo ""
    
    for T in 1 2 4 8; do
        echo "--- $T threads ---"
        ./build-ultra/bin/test_ultra_parallel $GRAPH $R $S $T 2>&1 | \
            grep -E "(Reference|Ultra Parallel|Speedup|Total time)"
        echo ""
    done
done

echo "=========================================="
echo "Server testing complete!"
echo "=========================================="
