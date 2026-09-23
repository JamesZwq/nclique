#!/bin/bash
# Build-time ablation on tods2 (2026-09-23): chain_index_tool --build-only under the four solver x tree-pass settings on
# the seven graphs of tods2.json and the eight vertex-induced samples of the scalability figure (scale/), two rounds.
set -u
cd ~/UNSW/pivoter
R=research/r1_skyline_index_20260918
G=/data/wenqianz/dsets
S=$R/scale
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg) git=$(git rev-parse --short HEAD)"
cmake --build $R/build -j 12 --target chain_index_tool || { echo "BUILD FAILED"; exit 1; }
python3 $R/run_build_ablation.py tods2_ablation 2 $G/com-amazon.edges $G/com-dblp.edges $G/web-NotreDame.edges $G/web-Stanford.edges $G/web-Google.edges \
  $S/cit-Patents_p20.edges $S/cit-Patents_p40.edges $S/cit-Patents_p60.edges $S/cit-Patents_p80.edges $G/cit-Patents.edges \
  $S/web-BerkStan_p20.edges $S/web-BerkStan_p40.edges $S/web-BerkStan_p60.edges $S/web-BerkStan_p80.edges $G/web-BerkStan.edges
echo "END $(date -Is) rc=$?"
