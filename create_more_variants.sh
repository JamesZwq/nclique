#!/bin/bash

# 创建更多优化变体

BASE_DIR="/Users/zhangwenqian/UNSW/pivoter"
VARIANTS_DIR="$BASE_DIR/optimization_variants"
SOURCE_FILE="$BASE_DIR/src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp"

echo "Creating additional optimization variants..."

# 变体 5: 减少 heap 更新频率（批量更新）
echo "Creating variant 5: Batch heap updates"
cat > $VARIANTS_DIR/variant_05_batch_heap_update.cpp << 'VARIANT5'
// 在这里插入批量 heap 更新的代码
// 暂时复制原文件
VARIANT5
cp $SOURCE_FILE $VARIANTS_DIR/variant_05_batch_heap_update.cpp

# 在 heap 更新部分添加批量处理
# 这需要手动编辑，先标记位置

# 变体 6: 优化 Phase A 的 schedule
echo "Creating variant 6: Optimize Phase A schedule"
cp $SOURCE_FILE $VARIANTS_DIR/variant_06_phase_a_schedule.cpp
# 修改 Phase A 的 schedule 参数
sed -i '' 's/#pragma omp for schedule(dynamic, 1)/#pragma omp for schedule(dynamic, 8)/g' \
    $VARIANTS_DIR/variant_06_phase_a_schedule.cpp

# 变体 7: 同时优化两个 schedule
echo "Creating variant 7: Optimize both schedules"
cp $SOURCE_FILE $VARIANTS_DIR/variant_07_both_schedules.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 1)/#pragma omp for schedule(dynamic, 8)/g' \
    $VARIANTS_DIR/variant_07_both_schedules.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(dynamic, 32)/g' \
    $VARIANTS_DIR/variant_07_both_schedules.cpp

# 变体 8: 使用 guided schedule for Phase B
echo "Creating variant 8: Guided schedule for Phase B"
cp $SOURCE_FILE $VARIANTS_DIR/variant_08_guided_phase_b.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(guided)/g' \
    $VARIANTS_DIR/variant_08_guided_phase_b.cpp

echo "Created additional variants!"
ls -lh $VARIANTS_DIR/ | tail -5
