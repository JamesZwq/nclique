#!/bin/bash
# tods1, queued behind the CND sweep (tmux session prior1): commit and push its records, then the S trees latency baseline
# (run_strees.py) on the twenty tods1 graphs, then commit and push that too.  Timing runs never overlap: this waits.
set -u
cd /data/wenqianz/pivoter_repo
export OMP_NUM_THREADS=1
G=/data/wenqianz/graphs
GIT="git -c user.name=Wenqian\ Zhang -c user.email=zhangwenqian6915@gmail.com"
while tmux has-session -t prior1 2>/dev/null; do sleep 120; done
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
git add -f src-r1index/scripts/prior_original_*.json src-r1index/scripts/tods1_prior.log
$GIT commit -qm "tods1: CND (original single-size) sweep records for the twenty tods1 graphs" || true
git pull -q --rebase origin main && git push -q origin main && echo "pushed prior records"
python3 research/r1_skyline_index_20260918/run_strees.py --tag tods1 \
  $G/ca-GrQc.edges $G/ca-CondMat.edges $G/ca-HepPh.edges $G/ca-AstroPh.edges $G/email-Eu-core.edges $G/dblp-core30.edges \
  $G/soc-Epinions1.edges $G/com-amazon.ungraph.edges $G/ca-MathSciNet.edges $G/com-dblp.edges $G/ca-dblp-2012.edges \
  $G/web-Stanford.edges $G/web-Google.edges $G/com-youtube.edges $G/web-uk-2005.edges $G/web-it-2004.edges \
  $G/ca-coauthors-dblp.edges $G/soc-pokec-relationships.edges $G/tech-as-skitter.edges $G/wiki-Talk.edges
echo "strees rc=$?"
git add -f research/r1_skyline_index_20260918/stages/index_vertices_tods1.json research/r1_skyline_index_20260918/stages/index-logs_vertices_tods1
$GIT commit -qm "tods1: S trees latency baseline (stage-2 vertices mode) on the twenty tods1 graphs" || true
git pull -q --rebase origin main && git push -q origin main && echo "pushed strees records"
echo "END $(date -Is)"
