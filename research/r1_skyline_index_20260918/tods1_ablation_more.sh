#!/bin/bash
# tods1 (2026-09-23): four more rounds (rounds 2-5) of the build-time ablation for the fourteen graphs below ten seconds;
# two rounds on this shared server left +-15% noise (web-Google final setting 5.21 s and 6.97 s).
set -u
cd /data/wenqianz/pivoter_repo
R=research/r1_skyline_index_20260918
G=/data/wenqianz/graphs
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg) git=$(git rev-parse --short HEAD)"
python3 $R/run_build_ablation.py tods1_ablation 6 $G/ca-GrQc.edges $G/dblp-core30.edges $G/ca-CondMat.edges $G/email-Eu-core.edges \
  $G/ca-AstroPh.edges $G/ca-HepPh.edges $G/ca-MathSciNet.edges $G/com-amazon.ungraph.edges $G/ca-dblp-2012.edges $G/com-dblp.edges \
  $G/web-Stanford.edges $G/web-Google.edges $G/com-youtube.edges $G/soc-Epinions1.edges
echo "END $(date -Is) rc=$?"
