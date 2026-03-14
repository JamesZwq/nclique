#!/bin/bash
export OMP_NUM_THREADS=16
WORK_DIR="/home/wenqianz/nclique"
RESULTS_FILE="$WORK_DIR/results_16threads.txt"
cd $WORK_DIR
echo "Optimization Tests - 16 Threads - $(date)" > $RESULTS_FILE
for variant in optimization_variants/variant_*.cpp; do
    name=$(basename $variant .cpp)
    echo "Testing $name..." | tee -a $RESULTS_FILE
    cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak
    cp $variant src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
    cd build && make -j8 degeneracy_cliques > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "  COMPILE FAILED" | tee -a $RESULTS_FILE
        cd $WORK_DIR
        cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
        continue
    fi
    cd $WORK_DIR
    for i in 1 2 3; do
        echo "  Run $i:" | tee -a $RESULTS_FILE
        timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep -E "(Multi-threading|speedup)" | tee -a $RESULTS_FILE
        sleep 2
    done
    echo "" | tee -a $RESULTS_FILE
    cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
done
echo "All tests completed - $(date)" >> $RESULTS_FILE
