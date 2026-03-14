#!/bin/bash

# 自动化优化实验框架
# 会在服务器上运行一整晚，测试所有优化方案

RESULTS_DIR="/home/wenqianz/nclique/optimization_results"
mkdir -p $RESULTS_DIR

# 记录函数
log_result() {
    local variant=$1
    local dataset=$2
    local time=$3
    echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant,$dataset,$time" >> $RESULTS_DIR/results.csv
}

# 测试函数
test_variant() {
    local variant_name=$1
    local variant_file=$2
    
    echo "========================================"
    echo "Testing: $variant_name"
    echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "========================================"
    
    cd /home/wenqianz/nclique
    
    # 复制变体文件
    cp "$variant_file" src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    
    # 编译
    cd build
    make -j8 degeneracy_cliques > compile.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "编译失败！查看 compile.log"
        log_result "$variant_name" "compile" "FAILED"
        cd ..
        return 1
    fi
    
    echo "编译成功！"
    
    # 测试 com-dblp (3次取平均)
    echo "Testing com-dblp (3 runs)..."
    local total_time=0
    local success_count=0
    
    for i in 1 2 3; do
        echo "  Run $i/3..."
        local output=$(timeout 180 ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1)
        local time=$(echo "$output" | grep "Multi-threading" | grep -oP '\d+\.\d+(?= ms)')
        
        if [ -n "$time" ]; then
            total_time=$(echo "$total_time + $time" | bc)
            success_count=$((success_count + 1))
            echo "    Time: $time ms"
        else
            echo "    Failed or timeout"
        fi
    done
    
    if [ $success_count -gt 0 ]; then
        local avg_time=$(echo "scale=2; $total_time / $success_count" | bc)
        echo "  Average: $avg_time ms"
        log_result "$variant_name" "com-dblp" "$avg_time"
    else
        log_result "$variant_name" "com-dblp" "FAILED"
    fi
    
    # 测试 web-Google (1次)
    echo "Testing web-Google..."
    local output=$(timeout 600 ./bin/degeneracy_cliques /data/wenqianz/web-Google.edges 3 4 2>&1)
    local time=$(echo "$output" | grep "Multi-threading" | grep -oP '\d+\.\d+(?= ms)')
    
    if [ -n "$time" ]; then
        echo "  Time: $time ms"
        log_result "$variant_name" "web-Google" "$time"
    else
        echo "  Failed or timeout"
        log_result "$variant_name" "web-Google" "FAILED"
    fi
    
    cd ..
    echo ""
}

# 主程序
echo "Starting optimization experiments..."
echo "Results will be saved to: $RESULTS_DIR/results.csv"
echo "variant,dataset,time" > $RESULTS_DIR/results.csv

# 测试基准版本
test_variant "baseline" "/home/wenqianz/nclique/src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp"

# 这里会添加更多变体...

echo "========================================"
echo "All experiments completed!"
echo "========================================"
echo "Results summary:"
cat $RESULTS_DIR/results.csv
