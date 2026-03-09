#!/bin/bash

# 测试不同的 schedule 参数
# 这个脚本会创建多个版本并测试

VARIANTS=(
    "dynamic,1:当前版本"
    "dynamic,4:dynamic-4"
    "dynamic,8:dynamic-8"
    "dynamic,16:dynamic-16"
    "dynamic,32:dynamic-32"
    "guided:guided"
    "static:static"
)

cd /home/wenqianz/nclique

for variant in "${VARIANTS[@]}"; do
    IFS=':' read -r schedule desc <<< "$variant"
    
    echo "========================================"
    echo "Testing: $desc (schedule=$schedule)"
    echo "========================================"
    
    # 修改代码
    sed -i "s/#pragma omp for schedule(dynamic, 16)/#pragma omp for schedule($schedule)/" \
        src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    
    # 编译
    cd build
    make -j8 degeneracy_cliques > /dev/null 2>&1
    
    if [ $? -ne 0 ]; then
        echo "编译失败！"
        continue
    fi
    
    # 测试 com-dblp
    echo "Testing com-dblp..."
    timeout 120 ./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | \
        grep -E "(Multi-threading|took:)" | tail -5
    
    echo ""
    cd ..
done

echo "========================================"
echo "所有测试完成！"
echo "========================================"
