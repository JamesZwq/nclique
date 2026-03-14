#!/bin/bash
export OMP_NUM_THREADS=16
cd /home/wenqianz/nclique
echo "Start: $(date)" > results_16t.txt
for v in optimization_variants/variant_*.cpp; do
  n=$(basename $v .cpp)
  echo "Test $n" | tee -a results_16t.txt
  cp $v src/NucleusDecomposition/NucleusCoreDecompositionRemoveScliqueRef.cpp
  cd build && make -j8 degeneracy_cliques >/dev/null 2>&1 && cd ..
  timeout 180 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep "Multi-threading" | tee -a results_16t.txt
  sleep 1
done
echo "End: $(date)" >> results_16t.txt
