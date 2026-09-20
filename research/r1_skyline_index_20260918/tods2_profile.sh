#!/bin/bash
# Query-latency profiles on tods2 (after the scalability session scale2 has finished): the seven graphs and the eight samples.
set -u
cd ~/UNSW/pivoter
R=research/r1_skyline_index_20260918
export OMP_NUM_THREADS=1
while tmux has-session -t scale2 2>/dev/null; do sleep 60; done
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
python3 $R/run_profile.py --tag profile_tods2 $R/cx/com-dblp.cx $R/cx/web-Stanford.cx $R/cx/com-amazon.cx $R/cx/web-NotreDame.cx $R/cx/web-Google.cx $R/cx/web-BerkStan.cx $R/cx/cit-Patents.cx \
  $R/cx/cit-Patents_p20.cx $R/cx/cit-Patents_p40.cx $R/cx/cit-Patents_p60.cx $R/cx/cit-Patents_p80.cx $R/cx/web-BerkStan_p20.cx $R/cx/web-BerkStan_p40.cx $R/cx/web-BerkStan_p60.cx $R/cx/web-BerkStan_p80.cx
echo "END $(date -Is) rc=$?"
