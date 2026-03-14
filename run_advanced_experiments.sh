#!/bin/bash

# 在服务器上运行完整的实验
# 测试1, 2, 4, 8线程的性能

set -e

echo "=========================================="
echo "Advanced Parallel Algorithm Experiments"
echo "Target: 3x+ speedup on 8 threads"
echo "=========================================="
echo ""

# 检测CPU信息
if [ -f /proc/cpuinfo ]; then
    CORES=$(grep -c ^processor /proc/cpuinfo)
    echo "CPU cores available: $CORES"
    grep "model name" /proc/cpuinfo | head -1
    echo ""
fi

# 测试数据集
DATASETS=(
    "toy.edges 3 4"
    "new_small_garph.edges 3 4"
)

# 如果有更大的数据集，添加到这里
if [ -f "data/com-dblp.ungraph.txt" ]; then
    DATASETS+=("data/com-dblp.ungraph.txt 3 4")
fi

if [ -f "data/web-Google.txt" ]; then
    DATASETS+=("data/web-Google.txt 3 4")
fi

# 线程数配置
THREAD_COUNTS=(1 2 4 8)

# 结果文件
RESULT_FILE="advanced_parallel_results_$(date +%Y%m%d_%H%M%S).txt"

echo "Results will be saved to: $RESULT_FILE"
echo ""

# 确保编译最新版本
echo "=========================================="
echo "Compiling..."
echo "=========================================="

if [ ! -f "build-advanced/test_advanced_parallel" ]; then
    ./build_and_test_advanced.sh
fi

echo ""
echo "=========================================="
echo "Starting Experiments"
echo "=========================================="
echo "" | tee -a "$RESULT_FILE"

# 对每个数据集运行实验
for dataset_config in "${DATASETS[@]}"; do
    read -r dataset r s <<< "$dataset_config"
    
    if [ ! -f "$dataset" ]; then
        echo "Skipping $dataset (file not found)"
        continue
    fi
    
    echo "========================================" | tee -a "$RESULT_FILE"
    echo "Dataset: $dataset" | tee -a "$RESULT_FILE"
    echo "Parameters: r=$r, s=$s" | tee -a "$RESULT_FILE"
    echo "========================================" | tee -a "$RESULT_FILE"
    echo "" | tee -a "$RESULT_FILE"
    
    # 对每个线程数运行测试
    for threads in "${THREAD_COUNTS[@]}"; do
        echo "----------------------------------------" | tee -a "$RESULT_FILE"
        echo "Testing with $threads thread(s)..." | tee -a "$RESULT_FILE"
        echo "----------------------------------------" | tee -a "$RESULT_FILE"
        
        # 运行3次取平均
        total_time=0
        total_speedup=0
        success_count=0
        
        for run in 1 2 3; do
            echo "Run $run/3..." | tee -a "$RESULT_FILE"
            
            # 运行测试
            output=$(./build-advanced/test_advanced_parallel "$dataset" "$r" "$s" "$threads" 2>&1)
            
            # 检查是否成功
            if echo "$output" | grep -q "ALL TESTS PASSED"; then
                # 提取时间和加速比
                ref_time=$(echo "$output" | grep "Reference (1 thread):" | awk '{print $4}')
                adv_time=$(echo "$output" | grep "Advanced ($threads threads):" | awk '{print $4}')
                speedup=$(echo "$output" | grep "Speedup:" | awk '{print $2}' | sed 's/x//')
                
                echo "  Reference: ${ref_time}ms, Advanced: ${adv_time}ms, Speedup: ${speedup}x" | tee -a "$RESULT_FILE"
                
                total_time=$(echo "$total_time + $adv_time" | bc)
                total_speedup=$(echo "$total_speedup + $speedup" | bc)
                success_count=$((success_count + 1))
            else
                echo "  ERROR: Test failed!" | tee -a "$RESULT_FILE"
                echo "$output" | tail -20 | tee -a "$RESULT_FILE"
            fi
        done
        
        # 计算平均值
        if [ $success_count -gt 0 ]; then
            avg_time=$(echo "scale=2; $total_time / $success_count" | bc)
            avg_speedup=$(echo "scale=3; $total_speedup / $success_count" | bc)
            
            echo "" | tee -a "$RESULT_FILE"
            echo "Average ($threads threads): ${avg_time}ms, Speedup: ${avg_speedup}x" | tee -a "$RESULT_FILE"
            
            # 检查是否达到目标
            if [ "$threads" -eq 8 ]; then
                target_met=$(echo "$avg_speedup >= 3.0" | bc)
                if [ "$target_met" -eq 1 ]; then
                    echo "✓ TARGET ACHIEVED: 3x+ speedup on 8 threads!" | tee -a "$RESULT_FILE"
                else
                    echo "✗ Target not met (need 3x+, got ${avg_speedup}x)" | tee -a "$RESULT_FILE"
                fi
            fi
        else
            echo "All runs failed for $threads threads" | tee -a "$RESULT_FILE"
        fi
        
        echo "" | tee -a "$RESULT_FILE"
    done
    
    echo "" | tee -a "$RESULT_FILE"
done

echo "========================================" | tee -a "$RESULT_FILE"
echo "Experiments Complete!" | tee -a "$RESULT_FILE"
echo "Results saved to: $RESULT_FILE" | tee -a "$RESULT_FILE"
echo "========================================" | tee -a "$RESULT_FILE"
