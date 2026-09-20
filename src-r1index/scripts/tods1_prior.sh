#!/bin/bash
# Original single-size implementation (NCliqueVertexCoreDecomposition) run once per size on every tods1 graph, small to large.
set -u
cd /data/wenqianz/pivoter_repo
export OMP_NUM_THREADS=1
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg)"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-GrQc.edges 44 ca-GrQc && echo "DONE ca-GrQc"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-CondMat.edges 26 ca-CondMat && echo "DONE ca-CondMat"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-HepPh.edges 239 ca-HepPh && echo "DONE ca-HepPh"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-AstroPh.edges 57 ca-AstroPh && echo "DONE ca-AstroPh"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/email-Eu-core.edges 35 email-Eu-core && echo "DONE email-Eu-core"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/dblp-core30.edges 114 dblp-core30 && echo "DONE dblp-core30"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/soc-Epinions1.edges 68 soc-Epinions1 && echo "DONE soc-Epinions1"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/com-amazon.ungraph.edges 7 com-amazon.ungraph && echo "DONE com-amazon.ungraph"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-MathSciNet.edges 25 ca-MathSciNet && echo "DONE ca-MathSciNet"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/com-dblp.edges 114 com-dblp && echo "DONE com-dblp"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-dblp-2012.edges 114 ca-dblp-2012 && echo "DONE ca-dblp-2012"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/web-Stanford.edges 72 web-Stanford && echo "DONE web-Stanford"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/web-Google.edges 45 web-Google && echo "DONE web-Google"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/com-youtube.edges 52 com-youtube && echo "DONE com-youtube"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/web-uk-2005.edges 500 web-uk-2005 && echo "DONE web-uk-2005"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/web-it-2004.edges 432 web-it-2004 && echo "DONE web-it-2004"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/ca-coauthors-dblp.edges 337 ca-coauthors-dblp && echo "DONE ca-coauthors-dblp"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/soc-pokec-relationships.edges 48 soc-pokec-relationships && echo "DONE soc-pokec-relationships"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/tech-as-skitter.edges 112 tech-as-skitter && echo "DONE tech-as-skitter"
python3 src-r1index/scripts/prior_sweep.py /data/wenqianz/graphs/wiki-Talk.edges 132 wiki-Talk && echo "DONE wiki-Talk"
echo "END $(date -Is)"
