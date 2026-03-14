#!/bin/bash

# 创建所有优化变体

BASE_DIR="/Users/zhangwenqian/UNSW/pivoter"
VARIANTS_DIR="$BASE_DIR/optimization_variants"
SOURCE_FILE="$BASE_DIR/src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp"

mkdir -p $VARIANTS_DIR

echo "Creating optimization variants..."

# 变体 1: schedule(dynamic, 4)
echo "Creating variant 1: schedule(dynamic, 4)"
cp $SOURCE_FILE $VARIANTS_DIR/variant_01_schedule_dynamic_4.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(dynamic, 4)/g' \
    $VARIANTS_DIR/variant_01_schedule_dynamic_4.cpp

# 变体 2: schedule(dynamic, 32)
echo "Creating variant 2: schedule(dynamic, 32)"
cp $SOURCE_FILE $VARIANTS_DIR/variant_02_schedule_dynamic_32.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(dynamic, 32)/g' \
    $VARIANTS_DIR/variant_02_schedule_dynamic_32.cpp

# 变体 3: schedule(guided)
echo "Creating variant 3: schedule(guided)"
cp $SOURCE_FILE $VARIANTS_DIR/variant_03_schedule_guided.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(guided)/g' \
    $VARIANTS_DIR/variant_03_schedule_guided.cpp

# 变体 4: schedule(static)
echo "Creating variant 4: schedule(static)"
cp $SOURCE_FILE $VARIANTS_DIR/variant_04_schedule_static.cpp
sed -i '' 's/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule(static)/g' \
    $VARIANTS_DIR/variant_04_schedule_static.cpp

echo "Created 4 variants in $VARIANTS_DIR"
ls -lh $VARIANTS_DIR/
