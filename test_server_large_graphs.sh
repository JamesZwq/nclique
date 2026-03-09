#!/bin/bash

# 服务器性能测试脚本 - 使用大图

echo "=========================================="
echo "Server Performance Testing (LARGE GRAPHS)"
echo "Multi-threaded Nucleus Decomposition"
echo "=========================================="
echo ""

# 检查环境
echo "Environment Information:"
echo "  Hostname: $(hostname)"
echo "  CPU Cores: $(nproc)"
echo "  Memory: $(free -h | grep Mem | awk '{print $2}')"
echo ""

# 设置测试配置
GRAPHS=("test_server_large_2000.edges" "test_server_xlarge_5000.edges")
GRAPH_NAMES=("Large (2000 nodes)" "Extra Large (5000 nodes)")
THREADS=(1 2 4 8)
R=3
S=4

# 结果文件
RESULTS_FILE="server_large_results_$(date +%Y%m%d_%H%M%S).txt"

echo "Test Configuration:" | tee $RESULTS_FILE
echo "  Graphs: ${GRAPH_NAMES[@]}" | tee -a $RESULTS_FILE
echo "  Thread counts: ${THREADS[@]}" | tee -a $RESULTS_FILE
echo "  Parameters: r=$R, s=$S" | tee -a $RESULTS_FILE
echo "" | tee -a $RESULTS_FILE
echo "========================================" | tee -a $RESULTS_FILE
echo "" | tee -a $RESULTS_FILE

# 测试每个图
for i in "${!GRAPHS[@]}"; do
    GRAPH="${GRAPHS[$i]}"
    GRAPH_NAME="${GRAPH_NAMES[$i]}"
    
    if [ ! -f "$GRAPH" ]; then
        echo "Warning: $GRAPH not found, skipping..." | tee -a $RESULTS_FILE
        continue
    fi
    
    echo "" | tee -a $RESULTS_FILE
    echo "========================================" | tee -a $RESULTS_FILE
    echo "Testing: $GRAPH_NAME" | tee -a $RESULTS_FILE
    echo "File: $GRAPH" | tee -a $RESULTS_FILE
    echo "========================================" | tee -a $RESULTS_FILE
    echo "" | tee -a $RESULTS_FILE
    
    # 存储基准时间
    BASELINE_TIME=""
    REF_TIME=""
    
    # 测试不同线程数
    for T in "${THREADS[@]}"; do
        echo "Testing with $T thread(s)..." | tee -a $RESULTS_FILE
        
        # 运行测试
        OUTPUT=$(./build-ultra/bin/test_ultra_parallel $GRAPH $R $S $T 2>&1)
        
        # 提取时间
        ULTRA_TIME=$(echo "$OUTPUT" | grep "Ultra Parallel ($T" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
        
        if [ $T -eq 1 ]; then
            REF=$(echo "$OUTPUT" | grep "Reference (1 thread)" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
            REF_TIME=$REF
            BASELINE_TIME=$ULTRA_TIME
            echo "  Reference (1 thread): ${REF}ms" | tee -a $RESULTS_FILE
            echo "  Ultra Parallel (1 thread): ${ULTRA_TIME}ms (baseline)" | tee -a $RESULTS_FILE
        else
            echo "  Ultra Parallel ($T threads): ${ULTRA_TIME}ms" | tee -a $RESULTS_FILE
            
            # 计算加速比
            if [ ! -z "$BASELINE_TIME" ] && [ $BASELINE_TIME -gt 0 ]; then
                SPEEDUP=$(echo "scale=2; $BASELINE_TIME / $ULTRA_TIME" | bc)
                echo "    Speedup vs 1 thread: ${SPEEDUP}x" | tee -a $RESULTS_FILE
            fi
            
            if [ ! -z "$REF_TIME" ] && [ $REF_TIME -gt 0 ]; then
                SPEEDUP_REF=$(echo "scale=2; $REF_TIME / $ULTRA_TIME" | bc)
                echo "    Speedup vs Reference: ${SPEEDUP_REF}x" | tee -a $RESULTS_FILE
            fi
        fi
        
        echo "" | tee -a $RESULTS_FILE
    done
done

echo "" | tee -a $RESULTS_FILE
echo "========================================" | tee -a $RESULTS_FILE
echo "Testing Complete!" | tee -a $RESULTS_FILE
echo "========================================" | tee -a $RESULTS_FILE
echo "" | tee -a $RESULTS_FILE
echo "Results saved to: $RESULTS_FILE" | tee -a $RESULTS_FILE
echo "" | tee -a $RESULTS_FILE

# 显示最佳结果
echo "Summary of Best Results:" | tee -a $RESULTS_FILE
echo "------------------------" | tee -a $RESULTS_FILE
grep "Speedup vs Reference:" $RESULTS_FILE | sort -t: -k2 -rn | head -5 | tee -a $RESULTS_FILE

echo ""
echo "For detailed results, see: $RESULTS_FILE"
