#!/bin/bash

# 系统性能测试脚本

set -e

echo "=========================================="
echo "Multi-threaded Performance Testing"
echo "=========================================="
echo ""

# 测试配置
GRAPHS=("test_small.edges" "test_medium.edges" "test_large.edges")
GRAPH_NAMES=("Small (50 nodes)" "Medium (200 nodes)" "Large (500 nodes)")
THREADS=(1 2 4 8)
R=3
S=4

# 结果文件
RESULTS_FILE="performance_results.txt"
echo "Performance Test Results - $(date)" > $RESULTS_FILE
echo "========================================" >> $RESULTS_FILE
echo "" >> $RESULTS_FILE

# 测试每个图
for i in "${!GRAPHS[@]}"; do
    GRAPH="${GRAPHS[$i]}"
    GRAPH_NAME="${GRAPH_NAMES[$i]}"
    
    if [ ! -f "$GRAPH" ]; then
        echo "Warning: $GRAPH not found, skipping..."
        continue
    fi
    
    echo ""
    echo "=========================================="
    echo "Testing: $GRAPH_NAME"
    echo "=========================================="
    echo ""
    echo "Graph: $GRAPH_NAME" >> $RESULTS_FILE
    echo "File: $GRAPH" >> $RESULTS_FILE
    echo "Parameters: r=$R, s=$S" >> $RESULTS_FILE
    echo "" >> $RESULTS_FILE
    
    # 测试不同线程数
    for T in "${THREADS[@]}"; do
        echo "Testing with $T thread(s)..."
        
        # 运行测试并提取时间
        OUTPUT=$(./build-ultra/bin/test_ultra_parallel $GRAPH $R $S $T 2>&1)
        
        # 提取Ultra Parallel时间
        ULTRA_TIME=$(echo "$OUTPUT" | grep "Ultra Parallel" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
        
        # 提取Reference时间（只在1线程时）
        if [ $T -eq 1 ]; then
            REF_TIME=$(echo "$OUTPUT" | grep "Reference" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
            echo "  Reference (1 thread): ${REF_TIME} ms" | tee -a $RESULTS_FILE
        fi
        
        echo "  Ultra Parallel ($T threads): ${ULTRA_TIME} ms" | tee -a $RESULTS_FILE
        
        # 计算加速比（相对于1线程Ultra Parallel）
        if [ $T -eq 1 ]; then
            BASELINE_TIME=$ULTRA_TIME
        else
            if [ $BASELINE_TIME -gt 0 ]; then
                SPEEDUP=$(echo "scale=2; $BASELINE_TIME / $ULTRA_TIME" | bc)
                echo "    Speedup vs 1 thread: ${SPEEDUP}x" | tee -a $RESULTS_FILE
            fi
        fi
    done
    
    echo "" >> $RESULTS_FILE
    echo "----------------------------------------" >> $RESULTS_FILE
    echo "" >> $RESULTS_FILE
done

echo ""
echo "=========================================="
echo "All tests complete!"
echo "Results saved to: $RESULTS_FILE"
echo "=========================================="
