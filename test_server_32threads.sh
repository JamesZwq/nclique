#!/bin/bash

# 服务器性能测试脚本 - 支持最多32线程

echo "=========================================="
echo "Server Performance Testing"
echo "Multi-threaded Nucleus Decomposition"
echo "=========================================="
echo ""

# 检查环境
echo "Environment Information:"
echo "  Hostname: $(hostname)"
echo "  CPU Info: $(grep "model name" /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)"
echo "  CPU Cores: $(nproc)"
echo "  Memory: $(free -h | grep Mem | awk '{print $2}')"
echo ""

# 设置测试配置
GRAPHS=("test_small.edges" "test_medium_500.edges" "test_large_1000.edges")
GRAPH_NAMES=("Small (50 nodes)" "Medium (500 nodes)" "Large (1000 nodes)")
THREADS=(1 2 4 8 16 32)
R=3
S=4

# 结果文件
RESULTS_FILE="server_performance_results_$(date +%Y%m%d_%H%M%S).txt"

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
        # 检查线程数是否超过CPU核心数
        MAX_THREADS=$(nproc)
        if [ $T -gt $MAX_THREADS ]; then
            echo "Skipping $T threads (exceeds available cores: $MAX_THREADS)" | tee -a $RESULTS_FILE
            continue
        fi
        
        echo "Testing with $T thread(s)..." | tee -a $RESULTS_FILE
        
        # 运行3次取平均
        TOTAL_ULTRA=0
        TOTAL_REF=0
        RUNS=3
        
        for RUN in $(seq 1 $RUNS); do
            OUTPUT=$(./build/bin/test_ultra_parallel $GRAPH $R $S $T 2>&1)
            
            # 提取时间
            ULTRA_TIME=$(echo "$OUTPUT" | grep "Ultra Parallel ($T" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
            TOTAL_ULTRA=$((TOTAL_ULTRA + ULTRA_TIME))
            
            if [ $T -eq 1 ]; then
                REF=$(echo "$OUTPUT" | grep "Reference (1 thread)" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
                TOTAL_REF=$((TOTAL_REF + REF))
            fi
        done
        
        # 计算平均值
        AVG_ULTRA=$((TOTAL_ULTRA / RUNS))
        
        if [ $T -eq 1 ]; then
            AVG_REF=$((TOTAL_REF / RUNS))
            REF_TIME=$AVG_REF
            BASELINE_TIME=$AVG_ULTRA
            echo "  Reference (1 thread): ${AVG_REF}ms" | tee -a $RESULTS_FILE
            echo "  Ultra Parallel (1 thread): ${AVG_ULTRA}ms (baseline)" | tee -a $RESULTS_FILE
        else
            echo "  Ultra Parallel ($T threads): ${AVG_ULTRA}ms" | tee -a $RESULTS_FILE
            
            # 计算加速比
            if [ ! -z "$BASELINE_TIME" ] && [ $BASELINE_TIME -gt 0 ]; then
                SPEEDUP=$(echo "scale=2; $BASELINE_TIME / $AVG_ULTRA" | bc)
                echo "    Speedup vs 1 thread: ${SPEEDUP}x" | tee -a $RESULTS_FILE
            fi
            
            if [ ! -z "$REF_TIME" ] && [ $REF_TIME -gt 0 ]; then
                SPEEDUP_REF=$(echo "scale=2; $REF_TIME / $AVG_ULTRA" | bc)
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
