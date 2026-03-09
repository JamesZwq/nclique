#!/bin/bash

# 详细的性能测试脚本

set -e

echo "=========================================="
echo "Comprehensive Performance Testing"
echo "=========================================="
echo ""

# 测试配置
GRAPH="test_medium_500.edges"
R=3
S=4
THREADS=(1 2 4 6 8)

echo "Graph: $GRAPH"
echo "Parameters: r=$R, s=$S"
echo ""
echo "=========================================="
echo ""

# 存储结果
declare -a ULTRA_TIMES
declare -a REF_TIME

# 测试每个线程数（运行3次取平均）
for T in "${THREADS[@]}"; do
    echo "Testing with $T thread(s)..."
    
    TOTAL_ULTRA=0
    TOTAL_REF=0
    RUNS=3
    
    for RUN in $(seq 1 $RUNS); do
        echo -n "  Run $RUN/$RUNS... "
        
        OUTPUT=$(./build-ultra/bin/test_ultra_parallel $GRAPH $R $S $T 2>&1)
        
        # 提取Ultra Parallel时间
        ULTRA_TIME=$(echo "$OUTPUT" | grep "Ultra Parallel" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
        TOTAL_ULTRA=$((TOTAL_ULTRA + ULTRA_TIME))
        
        # 提取Reference时间
        if [ $T -eq 1 ]; then
            REF=$(echo "$OUTPUT" | grep "Reference" | grep -oE '[0-9]+ ms' | grep -oE '[0-9]+')
            TOTAL_REF=$((TOTAL_REF + REF))
        fi
        
        echo "${ULTRA_TIME}ms"
    done
    
    # 计算平均值
    AVG_ULTRA=$((TOTAL_ULTRA / RUNS))
    ULTRA_TIMES[$T]=$AVG_ULTRA
    
    if [ $T -eq 1 ]; then
        AVG_REF=$((TOTAL_REF / RUNS))
        REF_TIME=$AVG_REF
        echo "  Average Reference: ${AVG_REF}ms"
    fi
    
    echo "  Average Ultra Parallel: ${AVG_ULTRA}ms"
    
    # 计算加速比
    if [ $T -gt 1 ]; then
        SPEEDUP=$(echo "scale=2; ${ULTRA_TIMES[1]} / $AVG_ULTRA" | bc)
        echo "  Speedup vs 1 thread: ${SPEEDUP}x"
        
        # 计算相对于Reference的加速比
        if [ $REF_TIME -gt 0 ]; then
            SPEEDUP_REF=$(echo "scale=2; $REF_TIME / $AVG_ULTRA" | bc)
            echo "  Speedup vs Reference: ${SPEEDUP_REF}x"
        fi
    fi
    
    echo ""
done

echo "=========================================="
echo "Summary"
echo "=========================================="
echo ""
echo "Reference (1 thread): ${REF_TIME}ms"
echo ""
echo "Ultra Parallel Performance:"
for T in "${THREADS[@]}"; do
    AVG=${ULTRA_TIMES[$T]}
    if [ $T -eq 1 ]; then
        echo "  $T thread:  ${AVG}ms (baseline)"
    else
        SPEEDUP=$(echo "scale=2; ${ULTRA_TIMES[1]} / $AVG" | bc)
        SPEEDUP_REF=$(echo "scale=2; $REF_TIME / $AVG" | bc)
        echo "  $T threads: ${AVG}ms (${SPEEDUP}x vs 1 thread, ${SPEEDUP_REF}x vs reference)"
    fi
done

echo ""
echo "=========================================="

# 检查是否达到目标
SPEEDUP_8=$(echo "scale=2; ${ULTRA_TIMES[1]} / ${ULTRA_TIMES[8]}" | bc)
TARGET=3.0

if (( $(echo "$SPEEDUP_8 >= $TARGET" | bc -l) )); then
    echo "✓ SUCCESS: 8-thread speedup (${SPEEDUP_8}x) >= target (${TARGET}x)"
else
    echo "✗ FAILED: 8-thread speedup (${SPEEDUP_8}x) < target (${TARGET}x)"
fi

echo "=========================================="
