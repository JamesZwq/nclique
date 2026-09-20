#!/bin/bash
# tods1, queued behind after1 (CND records + S trees baseline): query_profile on the twenty tods1 index files, then commit and push.
set -u
cd /data/wenqianz/pivoter_repo
R=research/r1_skyline_index_20260918
export OMP_NUM_THREADS=1
GIT="git -c user.name=Wenqian\ Zhang -c user.email=zhangwenqian6915@gmail.com"
while tmux has-session -t after1 2>/dev/null; do sleep 120; done
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
git pull -q --rebase origin main
python3 $R/run_profile.py --tag profile_tods1 $(ls $R/cx/*.cx)
echo "profile rc=$?"
git add -f $R/profile_tods1.json $R/profile_tods1-logs
$GIT commit -qm "tods1: query latency profiles (index against S trees) on the twenty tods1 graphs" || true
git pull -q --rebase origin main && git push -q origin main && echo "pushed profiles"
echo "END $(date -Is)"
