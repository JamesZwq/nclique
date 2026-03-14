#!/bin/bash
VARIANT=$1
WORK_DIR="/home/wenqianz/nclique"
cd $WORK_DIR
cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak
cp optimization_variants/${VARIANT}.cpp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
cd build && make -j8 degeneracy_cliques 2>&1 | tail -3
cd $WORK_DIR && timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep -E "(Multi-threading|took:)"
cp src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp.bak src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
