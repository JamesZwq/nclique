#!/bin/bash

# 修复版本：确保从正确目录运行

RESULTS_DIR="/home/wenqianz/nclique/optimization_results"
VARIANTS_DIR="/home/wenqianz/nclique/optimization_variants"
WORK_DIR="/home/wenqianz/nclique"

mkdir -p $RESULTS_DIR

echo "timestamp,variant,dataset,time_ms,status" > $RESULTS_DIR/results.csv

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
    echo "Compiling..."
    make -j8 degeneracy_cliques > compile.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "编译失败！"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,compile,0,FAILED" >> $RESULTS_DIR/results.csv
        cd $WORK_DIR
        cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak \
           src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
        return 1
    fi
    
    echo "编译成功！"
    
    # 回到工作目录（重要！nCr.txt 在这里）
    cd $WORK_DIR
    
    # 测试 com-dblp (3次)
    echo "Testing com-dblp (3 runs)..."
    local times=()
    
    for i in 1 2 3; do
        echo "  Run $i/3..."
        local output=$(timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1)
        
        # 提取时间
        local nucleus_time=$(echo "$output" | grep "Multi-threading" | grep -oP '\d+\.\d+' | head -1)
        
        if [ -n "$nucleus_time" ]; then
            times+=($nucleus_time)
            echo "    Time: $nucleus_time ms"
        else
            echo "    Failed"
        fi
        
        sleep 2
    done
    
    if [ ${#times[@]} -gt 0 ]; then
        local sum=0
        for t in "${times[@]}"; do
            sum=$(echo "$sum + $t" | bc)
        done
        local avg=$(echo "scale=2; $sum / ${#times[@]}" | bc)
        echo "  Average: $avg ms"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,com-dblp,$avg,SUCCESS" >> $RESULTS_DIR/results.csv
    else
        echo "  All failed"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,com-dblp,0,FAILED" >> $RESULTS_DIR/results.csv
    fi
    
    # 测试 web-Google
    echo "Testing web-Google..."
    local output=$(timeout 600 ./build/bin/degeneracy_cliques /data/wenqianz/web-Google.edges 3 4 2>&1)
    local nucleus_time=$(echo "$output" | grep "Multi-threading" | grep -oP '\d+\.\d+' | head -1)
    
    if [ -n "$nucleus_time" ]; then
        echo "  Time: $nucleus_time ms"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,web-Google,$nucleus_time,SUCCESS" >> $RESULTS_DIR/results.csv
    else
        echo "  Failed"
        echo "$(date '+%Y-%m-%d %H:%M:%S'),$variant_name,web-Google,0,FAILED" >> $RESULTS_DIR/results.csv
    fi
    
    # 恢复
    cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak \
       src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    
    echo ""
}

# 主程序
echo "========================================" | tee $RESULTS_DIR/experiment.log
echo "Starting Optimization Experiments (Fixed)" | tee -a $RESULTS_DIR/experiment.log
echo "Start time: $(date)" | tee -a $RESULTS_DIR/experiment.log
echo "========================================" | tee -a $RESULTS_DIR/experiment.log

# 测试基准
test_variant "baseline" "$WORK_DIR/src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp" | tee -a $RESULTS_DIR/experiment.log

# 测试所有变体
for variant_file in $VARIANTS_DIR/variant_*.cpp; do
    variant_name=$(basename $variant_file .cpp)
    test_variant "$variant_name" "$variant_file" | tee -a $RESULTS_DIR/experiment.log
done

echo "========================================" | tee -a $RESULTS_DIR/experiment.log
echo "All experiments completed!" | tee -a $RESULTS_DIR/experiment.log
echo "End time: $(date)" | tee -a $RESULTS_DIR/experiment.log
echo "========================================" | tee -a $RESULTS_DIR/experiment.log
cat $RESULTS_DIR/results.csv | tee -a $RESULTS_DIR/experiment.log
