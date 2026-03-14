#!/bin/bash

# 监控实验进度

echo "========================================="
echo "Optimization Experiments Monitor"
echo "========================================="
echo ""

# 检查实验是否在运行
echo "Checking if experiments are running..."
ssh tods2 "ps aux | grep run_all_experiments | grep -v grep"
echo ""

# 显示最新的日志
echo "Latest log output:"
echo "========================================="
ssh tods2 "tail -50 /home/wenqianz/nclique/experiment_output.log 2>/dev/null || echo 'Log file not found yet'"
echo ""

# 显示结果摘要
echo "========================================="
echo "Results so far:"
echo "========================================="
ssh tods2 "cat /home/wenqianz/nclique/optimization_results/results.csv 2>/dev/null || echo 'Results file not found yet'"
echo ""

# 分析结果
echo "========================================="
echo "Performance Analysis:"
echo "========================================="
ssh tods2 "
if [ -f /home/wenqianz/nclique/optimization_results/results.csv ]; then
    echo 'Baseline performance:'
    grep 'baseline' /home/wenqianz/nclique/optimization_results/results.csv | grep 'SUCCESS'
    echo ''
    echo 'Best performers:'
    grep 'SUCCESS' /home/wenqianz/nclique/optimization_results/results.csv | \
        grep -v 'baseline' | sort -t',' -k3 -n | head -5
else
    echo 'No results yet'
fi
"
