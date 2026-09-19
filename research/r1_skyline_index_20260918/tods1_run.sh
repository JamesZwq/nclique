#!/bin/bash
# Chain index experiments on tods1: latencies (run_final) on every available graph, small to large.
set -u
cd /data/wenqianz/pivoter_repo
G=/data/wenqianz/graphs
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
python3 research/r1_skyline_index_20260918/run_final.py --tag tods1 --only \
  $G/ca-GrQc.edges $G/ca-HepPh.edges $G/com-dblp.edges $G/web-Stanford.edges \
  $G/ca-AstroPh.edges $G/ca-CondMat.edges $G/soc-Epinions1.edges $G/email-Eu-core.edges $G/dblp-core30.edges \
  $G/com-amazon.ungraph.edges $G/web-Google.edges $G/ca-dblp-2012.edges $G/ca-MathSciNet.edges \
  $G/com-youtube.edges $G/soc-pokec-relationships.edges $G/tech-as-skitter.edges $G/wiki-Talk.edges \
  $G/ca-coauthors-dblp.edges $G/web-uk-2005.edges $G/web-it-2004.edges \
  $G/com-lj.edges $G/ca-hollywood-2009.edges $G/com-orkut.edges
echo "END $(date -Is) rc=$?"
