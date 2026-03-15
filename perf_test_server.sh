#!/bin/bash

# SDCT 性能测试脚本
# 测试不同版本和不同线程数的性能

set -e

GRAPH="/data/wenqianz/com-dblp.edges"
PROJECT_DIR="/home/wenqianz/nclique_tmp"
BUILD_DIR="$PROJECT_DIR/build"
BIN_DIR="$BUILD_DIR/bin"

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}SDCT Performance Test Script${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 检查图文件
if [ ! -f "$GRAPH" ]; then
    echo -e "${RED}Error: Graph file not found at $GRAPH${NC}"
    exit 1
fi

echo -e "${GREEN}Graph file: $GRAPH${NC}"
echo -e "${GREEN}Project directory: $PROJECT_DIR${NC}"
echo ""

# 进入项目目录
cd "$PROJECT_DIR"

# 检查是否需要编译
if [ ! -f "$BIN_DIR/degeneracy_cliques" ]; then
    echo -e "${YELLOW}Binary not found, compiling...${NC}"
    if [ ! -d "$BUILD_DIR" ]; then
        mkdir -p "$BUILD_DIR"
        cd "$BUILD_DIR"
        cmake .. -DCMAKE_BUILD_TYPE=Release
    fi
    cd "$BUILD_DIR"
    make degeneracy_cliques -j16
    cd "$PROJECT_DIR"
fi

echo -e "${GREEN}Binary ready: $BIN_DIR/degeneracy_cliques${NC}"
echo ""

# 测试函数
run_test() {
    local threads=$1
    local version=$2
    
    export OMP_NUM_THREADS=$threads
    
    echo -e "${YELLOW}Testing $version with $threads threads...${NC}"
    
    # 运行程序并捕获输出
    local output=$("$BIN_DIR/degeneracy_cliques" "$GRAPH" 2 3 2>&1)
    
    # 从输出中提取时间信息
    local time_line=$(echo "$output" | grep "Tree Build" | head -1)
    
    if [ -z "$time_line" ]; then
        echo -e "${RED}Failed to extract timing information${NC}"
        return 1
    fi
    
    echo "$time_line"
    echo ""
}

# 测试不同线程数
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Performance Test Results${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 创建结果文件
RESULT_FILE="/tmp/sdct_perf_results.txt"
> "$RESULT_FILE"

echo "SDCT Performance Test Results" >> "$RESULT_FILE"
echo "Graph: $GRAPH" >> "$RESULT_FILE"
echo "Date: $(date)" >> "$RESULT_FILE"
echo "" >> "$RESULT_FILE"
echo "Threads | Time (ms) | Speedup" >> "$RESULT_FILE"
echo "--------|-----------|--------" >> "$RESULT_FILE"

# 测试单线程（基准）
echo -e "${BLUE}Running baseline (1 thread)...${NC}"
baseline_output=$(run_test 1 "SDCT_Par5")
baseline_time=$(echo "$baseline_output" | grep -oP '\d+(?= ms)' | head -1)

if [ -z "$baseline_time" ]; then
    echo -e "${RED}Failed to get baseline time${NC}"
    exit 1
fi

echo -e "${GREEN}Baseline time: ${baseline_time}ms${NC}"
echo "1 | $baseline_time | 1.00x" >> "$RESULT_FILE"
echo ""

# 测试不同线程数
for threads in 2 4 8 16 32 64; do
    echo -e "${BLUE}Testing with $threads threads...${NC}"
    
    output=$(run_test $threads "SDCT_Par5")
    time_ms=$(echo "$output" | grep -oP '\d+(?= ms)' | head -1)
    
    if [ -z "$time_ms" ]; then
        echo -e "${RED}Failed to get time for $threads threads${NC}"
        continue
    fi
    
    # 计算加速比
    speedup=$(echo "scale=2; $baseline_time / $time_ms" | bc)
    
    echo -e "${GREEN}Time: ${time_ms}ms, Speedup: ${speedup}x${NC}"
    echo "$threads | $time_ms | ${speedup}x" >> "$RESULT_FILE"
    echo ""
done

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Test Results Summary${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
cat "$RESULT_FILE"
echo ""

# 保存结果到文件
echo -e "${GREEN}Results saved to: $RESULT_FILE${NC}"

# 显示完整的结果文件
echo ""
echo -e "${BLUE}Full Results:${NC}"
cat "$RESULT_FILE"
