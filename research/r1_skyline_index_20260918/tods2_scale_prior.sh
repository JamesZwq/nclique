#!/bin/bash
# CND (original single-size implementation) once per size on the eight vertex-induced samples, after the profile session
# (profile2) has finished; s_max of each sample is read from scale_tods2.json.  Records: src-r1index/scripts/prior_original_<sample>.json.
set -u
cd ~/UNSW/pivoter
R=research/r1_skyline_index_20260918
export OMP_NUM_THREADS=1
while tmux has-session -t profile2 2>/dev/null; do sleep 60; done
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
for g in cit-Patents_p20 cit-Patents_p40 cit-Patents_p60 cit-Patents_p80 web-BerkStan_p20 web-BerkStan_p40 web-BerkStan_p60 web-BerkStan_p80; do
  smax=$(python3 -c "
import json,sys
for r in json.load(open('$R/scale_tods2.json'))['runs']:
    if r['graph'].endswith('/$g.edges') and 'result' in r: print(r['result']['s_max'])")
  [ -z "$smax" ] && { echo "no s_max for $g"; continue; }
  python3 src-r1index/scripts/prior_sweep.py $R/scale/$g.edges $smax $g && echo "DONE $g"
done
echo "END $(date -Is)"
