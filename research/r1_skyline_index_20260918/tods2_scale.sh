#!/bin/bash
# Scalability experiment on tods2: vertex-induced samples (20/40/60/80 percent) of cit-Patents and web-BerkStan, then the
# full chain index bench (run_final.py) on the samples.  Waits for the S trees baseline session (strees2) to finish first.
set -u
cd ~/UNSW/pivoter
R=research/r1_skyline_index_20260918
G=/data/wenqianz/dsets
S=$R/scale
export OMP_NUM_THREADS=1
while tmux has-session -t strees2 2>/dev/null; do sleep 60; done
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
mkdir -p $S
python3 $R/sample_vertices.py $G/cit-Patents.edges $S 0.2 0.4 0.6 0.8
python3 $R/sample_vertices.py $G/web-BerkStan.edges $S 0.2 0.4 0.6 0.8
python3 $R/run_final.py --tag scale_tods2 --only \
  $S/cit-Patents_p20.edges $S/cit-Patents_p40.edges $S/cit-Patents_p60.edges $S/cit-Patents_p80.edges \
  $S/web-BerkStan_p20.edges $S/web-BerkStan_p40.edges $S/web-BerkStan_p60.edges $S/web-BerkStan_p80.edges
echo "END $(date -Is) rc=$?"
