#!/bin/bash

# 持续优化实验 - 不断测试、分析、改进

RESULTS_DIR="/home/wenqianz/nclique/optimization_results"
VARIANTS_DIR="/home/wenqianz/nclique/optimization_variants"
WORK_DIR="/home/wenqianz/nclique"

mkdir -p $RESULTS_DIR

# 初始化
echo "timestamp,variant,dataset,time_ms,status" > $RESULTS_DIR/results.csv
echo "========================================" > $RESULTS_DIR/analysis.log
echo "Continuous Optimization Experiments" >> $RESULTS_DIR/analysis.log
echo "Start: $(date)" >> $RESULTS_DIR/analysis.log
echo "========================================" >> $RESULTS_DIR/analysis.log

# 测试单个变体
test_variant() {
    local variant_name=$1
    local variant_file=$2
    
    echo "========================================"
    echo "Testing: $variant_name"
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================"
    
    cd $WORK_DIR
    
    # 备份
    cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp \
       src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak
    
    # 复制变体
    cp "$variant_file" src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    
    # 编译
    cd build
    make -j8 degeneracy_cliques > compile.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "编译失败"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,compile,0,FAILED" >> $RESULTS_DIR/results.csv
        cd $WORK_DIR
        cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak \
           src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
        return 1
    fi
    
    cd $WORK_DIR
    
    # 测试 com-dblp (3次)
    echo "Testing com-dblp..."
    local times=()
    
    for i in 1 2 3; do
        local output=$(timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1)
        local time=$(echo "$output" | grep "Multi-threading" | grep -oP '\d+\.\d+' | head -1)
        
        if [ -n "$time" ]; then
            times+=($time)
            echo "  Run $i: $time ms"
        fi
        sleep 1
    done
    
    if [ ${#times[@]} -gt 0 ]; then
        local sum=0
        for t in "${times[@]}"; do
            sum=$(echo "$sum + $t" | bc)
        done
        local avg=$(echo "scale=2; $sum / ${#times[@]}" | bc)
        echo "  Average: $avg ms"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,com-dblp,$avg,SUCCESS" >> $RESULTS_DIR/results.csv
        
        # 如果性能好，测试 web-Google
        if (( $(echo "$avg < 450" | bc -l) )); then
            echo "Good performance! Testing web-Google..."
            local output=$(timeout 600 ./build/bin/degeneracy_cliques /data/wenqianz/web-Google.edges 3 4 2>&1)
            local time=$(echo "$output" | grep "Multi-threading" | grep -oP '\d+\.\d+' | head -1)
            
            if [ -n "$time" ]; then
                echo "  web-Google: $time ms"
                echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,web-Google,$time,SUCCESS" >> $RESULTS_DIR/results.csv
            fi
        fi
    else
        echo "  All failed"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,com-dblp,0,FAILED" >> $RESULTS_DIR/results.csv
    fi
    
    # 恢复
    cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak \
       src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    
    echo ""
}

# 分析结果
analyze_results() {
    echo "========================================"
    echo "Current Results Analysis"
    echo "Time: $(date)"
    echo "========================================"
    
    if [ ! -f $RESULTS_DIR/results.csv ]; then
        echo "No results yet"
        return
    fi
    
    # 找出最佳性能
    echo "Best performers (com-dblp):"
    grep "com-dblp.*SUCCESS" $RESULTS_DIR/results.csv | \
        sort -t',' -k3 -n | head -5 | \
        awk -F',' '{printf "  %s: %.2f ms\n", $2, $3}'
    
    echo ""
    echo "Best performers (web-Google):"
    grep "web-Google.*SUCCESS" $RESULTS_DIR/results.csv | \
        sort -t',' -k3 -n | head -5 | \
        awk -F',' '{printf "  %s: %.2f ms\n", $2, $3}'
    
    echo ""
    
    # 保存分析
    echo "--- Analysis at $(date) ---" >> $RESULTS_DIR/analysis.log
    grep "SUCCESS" $RESULTS_DIR/results.csv | tail -10 >> $RESULTS_DIR/analysis.log
    echo "" >> $RESULTS_DIR/analysis.log
}

# 主循环
echo "Starting continuous optimization..."

# 第一轮：测试所有变体
echo "=== Round 1: Testing all variants ==="
for variant_file in $VARIANTS_DIR/variant_*.cpp; do
    variant_name=$(basename $variant_file .cpp)
    test_variant "$variant_name" "$variant_file"
    
    # 每测试 5 个变体，分析一次
    if [ $(($(grep -c "SUCCESS" $RESULTS_DIR/results.csv 2>/dev/null || echo 0) % 5)) -eq 0 ]; then
        analyze_results
    fi
done

# 最终分析
echo "========================================"
echo "Final Analysis"
echo "========================================"
analyze_results

echo ""
echo "All experiments completed!"
echo "Results saved to: $RESULTS_DIR/results.csv"
