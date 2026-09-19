#!/bin/bash
# Chain index experiments on tods2 (quiet machine): latencies (run_final) on the graphs present in /data/wenqianz/dsets.
set -u
cd ~/UNSW/pivoter
G=/data/wenqianz/dsets
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
python3 research/r1_skyline_index_20260918/run_final.py --tag tods2 --only \
  $G/com-dblp.edges $G/web-Stanford.edges $G/com-amazon.edges $G/web-NotreDame.edges $G/web-Google.edges $G/web-BerkStan.edges $G/cit-Patents.edges
echo "END $(date -Is) rc=$?"
