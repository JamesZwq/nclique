#!/bin/bash

# 创建高级优化变体

BASE_DIR="/Users/zhangwenqian/UNSW/pivoter"
VARIANTS_DIR="$BASE_DIR/optimization_variants"
SOURCE_FILE="$BASE_DIR/src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp"

echo "Creating advanced optimization variants..."

# 变体 9: 减少 atomic 操作的内存顺序（使用 relaxed）
echo "Creating variant 9: Relaxed memory order for atomics"
cp $SOURCE_FILE $VARIANTS_DIR/variant_09_relaxed_atomics.cpp
# 已经使用 relaxed，这个变体保持不变作为对照

# 变体 10: 增加 Phase B 的并行粒度
echo "Creating variant 10: Larger chunk size for Phase B"
cp $SOURCE_FILE $VARIANTS_DIR/variant_10_larger_chunk.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(dynamic, 64)/g' \
    $VARIANTS_DIR/variant_10_larger_chunk.cpp

# 变体 11: 减小 Phase B 的并行粒度
echo "Creating variant 11: Smaller chunk size for Phase B"
cp $SOURCE_FILE $VARIANTS_DIR/variant_11_smaller_chunk.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(dynamic, 8)/g' \
    $VARIANTS_DIR/variant_11_smaller_chunk.cpp

# 变体 12: 使用 static schedule 但指定 chunk size
echo "Creating variant 12: Static schedule with chunk"
cp $SOURCE_FILE $VARIANTS_DIR/variant_12_static_chunk.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(static, 32)/g' \
    $VARIANTS_DIR/variant_12_static_chunk.cpp

# 变体 13: 优化 Phase A 使用更大的 chunk
echo "Creating variant 13: Larger chunk for Phase A"
cp $SOURCE_FILE $VARIANTS_DIR/variant_13_phase_a_large.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 1)/#pragma omp for schedule(dynamic, 4)/g' \
    $VARIANTS_DIR/variant_13_phase_a_large.cpp

# 变体 14: 两个阶段都使用 guided
echo "Creating variant 14: Both phases guided"
cp $SOURCE_FILE $VARIANTS_DIR/variant_14_both_guided.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 1)/#pragma omp for schedule(guided)/g' \
    $VARIANTS_DIR/variant_14_both_guided.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(guided)/g' \
    $VARIANTS_DIR/variant_14_both_guided.cpp

# 变体 15: Phase A guided, Phase B dynamic
echo "Creating variant 15: Mixed schedules"
cp $SOURCE_FILE $VARIANTS_DIR/variant_15_mixed_schedule.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 1)/#pragma omp for schedule(guided)/g' \
    $VARIANTS_DIR/variant_15_mixed_schedule.cpp

# 变体 16: 优化 treeGraphV 更新的 chunk size
echo "Creating variant 16: Optimize treeGraphV chunk"
cp $SOURCE_FILE $VARIANTS_DIR/variant_16_treegraphv_chunk.cpp
sed -i '' 's/#pragma omp parallel for schedule(dynamic, 1024)/#pragma omp parallel for schedule(dynamic, 512)/g' \
    $VARIANTS_DIR/variant_16_treegraphv_chunk.cpp

echo "Created 8 advanced variants!"
ls -lh $VARIANTS_DIR/variant_[0-9][0-9]*.cpp | wc -l
