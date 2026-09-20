#!/bin/bash
# Original single-size implementation run once per size on every tods2 graph.
set -u
cd ~/UNSW/pivoter
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/com-dblp.edges 114 com-dblp && echo "DONE com-dblp"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/web-Stanford.edges 72 web-Stanford && echo "DONE web-Stanford"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/com-amazon.edges 7 com-amazon && echo "DONE com-amazon"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/web-NotreDame.edges 156 web-NotreDame && echo "DONE web-NotreDame"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/web-Google.edges 45 web-Google && echo "DONE web-Google"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/web-BerkStan.edges 202 web-BerkStan && echo "DONE web-BerkStan"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/dsets/cit-Patents.edges 65 cit-Patents && echo "DONE cit-Patents"
echo "END $(date -Is)"
