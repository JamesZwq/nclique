#!/bin/bash
# Build-time ablation on tods1 (2026-09-23): chain_index_tool --build-only under the four solver x tree-pass settings on
# the twenty graphs of tods1.json, small to large; one round for all, a second round for the fourteen below one minute.
set -u
cd /data/wenqianz/pivoter_repo
R=research/r1_skyline_index_20260918
G=/data/wenqianz/graphs
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg) git=$(git rev-parse --short HEAD)"
cmake --build $R/build -j 12 --target chain_index_tool || { echo "BUILD FAILED"; exit 1; }
SMALL="$G/ca-GrQc.edges $G/dblp-core30.edges $G/ca-CondMat.edges $G/email-Eu-core.edges $G/ca-AstroPh.edges $G/ca-HepPh.edges \
  $G/ca-MathSciNet.edges $G/com-amazon.ungraph.edges $G/ca-dblp-2012.edges $G/com-dblp.edges $G/web-Stanford.edges $G/web-Google.edges \
  $G/com-youtube.edges $G/soc-Epinions1.edges"
LARGE="$G/web-it-2004.edges $G/web-uk-2005.edges $G/soc-pokec-relationships.edges $G/ca-coauthors-dblp.edges $G/wiki-Talk.edges $G/tech-as-skitter.edges"
python3 $R/run_build_ablation.py tods1_ablation 1 $SMALL $LARGE
python3 $R/run_build_ablation.py tods1_ablation 2 $SMALL
echo "END $(date -Is) rc=$?"
