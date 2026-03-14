#!/bin/bash
cd /home/wenqianz/nclique
echo 'SINGLE THREAD:'
export OMP_NUM_THREADS=1
./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep 'Multi-threading'
echo ''
echo '16 THREADS:'
export OMP_NUM_THREADS=16
./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep 'Multi-threading'
