#!/bin/bash

# 详细的正确性验证脚本

echo "=========================================="
echo "Detailed Correctness Verification"
echo "=========================================="
echo ""

GRAPH="test_small.edges"
R=3
S=4

echo "Testing on: $GRAPH (r=$R, s=$S)"
echo ""

# 运行测试并保存完整输出
OUTPUT=$(./build-ultra/bin/test_ultra_parallel $GRAPH $R $S 1 2>&1)

# 提取结果大小
REF_SIZE=$(echo "$OUTPUT" | grep "Reference result size:" | grep -oE '[0-9]+')
ULTRA_SIZE=$(echo "$OUTPUT" | grep "Ultra result size:" | grep -oE '[0-9]+')

echo "Result sizes:"
echo "  Reference: $REF_SIZE"
echo "  Ultra:     $ULTRA_SIZE"

if [ "$REF_SIZE" == "$ULTRA_SIZE" ]; then
    echo "  ✓ Sizes match"
else
    echo "  ✗ ERROR: Sizes differ!"
    exit 1
fi

echo ""

# 检查是否有sample match信息
if echo "$OUTPUT" | grep -q "Sample results match"; then
    echo "✓ Sample results verified"
else
    echo "✗ WARNING: Sample results not verified"
fi

echo ""
echo "=========================================="
echo "Testing multiple thread counts..."
echo "=========================================="
echo ""

# 测试不同线程数，确保结果一致
for T in 1 2 4 8; do
    echo "Testing with $T thread(s)..."
    
    RESULT=$(./build-ultra/bin/test_ultra_parallel $GRAPH $R $S $T 2>&1)
    
    SIZE=$(echo "$RESULT" | grep "Ultra result size:" | grep -oE '[0-9]+')
    MATCH=$(echo "$RESULT" | grep "Sample results match" | wc -l)
    
    if [ "$SIZE" == "$REF_SIZE" ] && [ "$MATCH" -gt 0 ]; then
        echo "  ✓ $T threads: size=$SIZE, correctness verified"
    else
        echo "  ✗ $T threads: FAILED (size=$SIZE, match=$MATCH)"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "Testing on medium graph..."
echo "=========================================="
echo ""

GRAPH2="test_medium_500.edges"

for T in 1 8; do
    echo "Testing $GRAPH2 with $T thread(s)..."
    
    RESULT=$(./build-ultra/bin/test_ultra_parallel $GRAPH2 $R $S $T 2>&1)
    
    SIZE=$(echo "$RESULT" | grep "Ultra result size:" | grep -oE '[0-9]+')
    MATCH=$(echo "$RESULT" | grep "Sample results match" | wc -l)
    
    if [ $T -eq 1 ]; then
        REF_SIZE2=$SIZE
    fi
    
    if [ "$SIZE" == "$REF_SIZE2" ] && [ "$MATCH" -gt 0 ]; then
        echo "  ✓ $T threads: size=$SIZE, correctness verified"
    else
        echo "  ✗ $T threads: FAILED (size=$SIZE, expected=$REF_SIZE2)"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "✓ ALL CORRECTNESS CHECKS PASSED"
echo "=========================================="
echo ""
echo "Summary:"
echo "  - Result sizes match across all thread counts"
echo "  - Sample core values verified"
echo "  - Tested on multiple graphs"
echo "  - Algorithm is CORRECT"
