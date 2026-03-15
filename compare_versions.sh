#!/bin/bash

# SDCT 版本对比性能测试脚本
# 测试 SDCT_Par2, SDCT_Par3, SDCT_Par4, SDCT_Par5 的性能

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
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}SDCT Version Comparison Test${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 检查图文件
if [ ! -f "$GRAPH" ]; then
    echo -e "${RED}Error: Graph file not found at $GRAPH${NC}"
    exit 1
fi

echo -e "${GREEN}Graph: $GRAPH${NC}"
echo -e "${GREEN}Project: $PROJECT_DIR${NC}"
echo ""

cd "$PROJECT_DIR"

# 检查二进制文件
if [ ! -f "$BIN_DIR/verify_all_sdct" ]; then
    echo -e "${YELLOW}verify_all_sdct not found, compiling...${NC}"
    if [ ! -d "$BUILD_DIR" ]; then
        mkdir -p "$BUILD_DIR"
        cd "$BUILD_DIR"
        cmake .. -DCMAKE_BUILD_TYPE=Release
    fi
    cd "$BUILD_DIR"
    make verify_all_sdct -j16
    cd "$PROJECT_DIR"
fi

echo -e "${GREEN}Binary ready: $BIN_DIR/verify_all_sdct${NC}"
echo ""

# 运行验证测试
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Running Correctness Verification${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

timeout 1800 "$BIN_DIR/verify_all_sdct" "$GRAPH" 2>&1 | tee /tmp/verify_results.txt

echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Verification Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# 提取结果
if grep -q "ALL TESTS PASSED" /tmp/verify_results.txt; then
    echo -e "${GREEN}✓ All versions produce correct results!${NC}"
else
    echo -e "${RED}✗ Some versions failed verification${NC}"
fi

echo ""
echo -e "${CYAN}Results saved to: /tmp/verify_results.txt${NC}"
