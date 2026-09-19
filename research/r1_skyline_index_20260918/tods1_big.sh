#!/bin/bash
# Second tods1 run after the 64-bit member offsets: the three graphs that overflowed the 32-bit clique-tree index.
set -u
cd /data/wenqianz/pivoter_repo
G=/data/wenqianz/graphs
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
python3 research/r1_skyline_index_20260918/run_final.py --tag tods1_big --only $G/com-lj.edges $G/ca-hollywood-2009.edges $G/com-orkut.edges
echo "END $(date -Is) rc=$?"
