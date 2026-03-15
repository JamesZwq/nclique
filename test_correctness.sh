#!/bin/bash

# 测试 SDCT 和 SDCT_Parallel 的正确性对比

GRAPH="graphs/com-dblp.edges"
BIN="./build/bin/degeneracy_cliques"

echo "=========================================="
echo "Correctness Test: SDCT vs SDCT_Parallel"
echo "=========================================="
echo ""

# 测试 SDCT（1线程）
echo "Testing SDCT (baseline, 1 thread)..."
export OMP_NUM_THREADS=1
RESULT_SDCT=$("$BIN" "$GRAPH" 2 3 2>&1 | grep "num Leaf" | tail -1)
echo "SDCT result: $RESULT_SDCT"
echo ""

# 测试 SDCT_Parallel（1线程）
echo "Testing SDCT_Parallel (1 thread)..."
export OMP_NUM_THREADS=1
RESULT_PAR1=$("$BIN" "$GRAPH" 2 3 2>&1 | grep "num Leaf" | tail -1)
echo "SDCT_Parallel (1T) result: $RESULT_PAR1"
echo ""

# 测试 SDCT_Parallel（4线程）
echo "Testing SDCT_Parallel (4 threads)..."
export OMP_NUM_THREADS=4
RESULT_PAR4=$("$BIN" "$GRAPH" 2 3 2>&1 | grep "num Leaf" | tail -1)
echo "SDCT_Parallel (4T) result: $RESULT_PAR4"
echo ""

# 测试 SDCT_Parallel（8线程）
echo "Testing SDCT_Parallel (8 threads)..."
export OMP_NUM_THREADS=8
RESULT_PAR8=$("$BIN" "$GRAPH" 2 3 2>&1 | grep "num Leaf" | tail -1)
echo "SDCT_Parallel (8T) result: $RESULT_PAR8"
echo ""

echo "=========================================="
echo "Verification"
echo "=========================================="

if [ "$RESULT_SDCT" = "$RESULT_PAR1" ]; then
    echo "✓ SDCT_Parallel (1T) matches SDCT"
else
    echo "✗ SDCT_Parallel (1T) DIFFERS from SDCT"
    echo "  SDCT:           $RESULT_SDCT"
    echo "  SDCT_Parallel:  $RESULT_PAR1"
fi

if [ "$RESULT_PAR1" = "$RESULT_PAR4" ]; then
    echo "✓ SDCT_Parallel (4T) matches (1T)"
else
    echo "✗ SDCT_Parallel (4T) DIFFERS from (1T)"
fi

if [ "$RESULT_PAR1" = "$RESULT_PAR8" ]; then
    echo "✓ SDCT_Parallel (8T) matches (1T)"
else
    echo "✗ SDCT_Parallel (8T) DIFFERS from (1T)"
fi
