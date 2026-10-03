#!/bin/bash
# tods1 (2026-10-04): the chain-index pipelines of the paper on three new graphs, strictly one timing process at a time:
# (1) run_final --bench (sizes, ratio, latencies)       -> tods1_new3.json, tods1_new3-logs/, cx/<g>.cx
# (2) omega (largest clique) read from each stored .cx  -> tods1_new3_omega.log
# (3) build-time ablation, four settings, rounds as tods1_ablation (6 below 10 s, 2 below 60 s, else 1)
#                                                       -> tods1_new3_ablation.json, tods1_new3_ablation-logs/
# (4) query_profile on the three .cx                    -> profile_tods1_new3.json, profile_tods1_new3-logs/
# (5) CND (original NCliqueVertexCoreDecomposition, default r=1 path, no PIVOTER_* env), sizes 2 .. omega
#                                                       -> src-r1index/scripts/prior_original_<g>.json (+ copy in prior/tods1/)
# Every pipeline under `timeout 6h`; a failure is recorded and the script moves on.
set -u
cd /data/wenqianz/pivoter_repo
R=research/r1_skyline_index_20260918
G=/data/wenqianz/graphs
GRAPHS="inf-roadNet-CA sc-msdoor web-baidu-baike"
export OMP_NUM_THREADS=1
T="timeout 21600"
echo "START $(date -Is) host=$(hostname) load=$(cut -d' ' -f1-3 /proc/loadavg) git=$(git rev-parse --short HEAD)"

echo "== (1) bench $(date -Is)"
$T python3 $R/run_final.py --tag tods1_new3 --only $G/inf-roadNet-CA.edges $G/sc-msdoor.edges $G/web-baidu-baike.edges
echo "bench rc=$? $(date -Is)"

echo "== (2) omega $(date -Is)"
OD=/data/wenqianz/omega_reader; mkdir -p $OD
cat > $OD/omega.cpp <<'CPP'
// largest chain omega of a stored ChainIndex (.cx), read with ChainIndex<double>::load, as in the 2026-09-24 omega.json audit
#include "/data/wenqianz/pivoter_repo/src-r1index/chain_index.hpp"
#include <iostream>
int main(int argc, char** argv) { using namespace chainindex;
    for (int i = 1; i < argc; ++i) { const auto ix = ChainIndex<double>::load(argv[i]); int w = 0;
        for (uint32_t c = 0; c < ix.chains; ++c) w = std::max<int>(w, ix.omega[c]);
        std::cout << argv[i] << " " << w << " max_size " << ix.max_size << "\n"; } }
CPP
g++ -O2 -std=c++23 -I/data/wenqianz/pivoter_repo/src-r1index -o $OD/omega $OD/omega.cpp 2>&1 | tail -5
for g in $GRAPHS; do [ -f $R/cx/$g.cx ] && $OD/omega $R/cx/$g.cx; done | tee $R/tods1_new3_omega.log

echo "== (3) ablation $(date -Is)"
for g in $GRAPHS; do
  ROUNDS=$(python3 - "$R/tods1_new3.json" "$g" <<'PY'
import json, sys
rec = json.load(open(sys.argv[1])); g = sys.argv[2]
r = next((r for r in rec['runs'] if r['graph'].endswith('/' + g + '.edges') and 'result' in r), None)
if r is None: print(1); raise SystemExit
x = r['result']; t = (x['ti_ms'] + x['build_ms'] + x['compact_ms']) / 1000
print(6 if t < 10 else 2 if t < 60 else 1)
PY
)
  echo "ablation $g rounds=$ROUNDS $(date -Is)"
  $T python3 $R/run_build_ablation.py tods1_new3_ablation $ROUNDS $G/$g.edges
  echo "ablation $g rc=$? $(date -Is)"
done

echo "== (4) profile $(date -Is)"
CX=""; for g in $GRAPHS; do [ -f $R/cx/$g.cx ] && CX="$CX $R/cx/$g.cx"; done
$T python3 $R/run_profile.py --tag profile_tods1_new3 $CX
echo "profile rc=$? $(date -Is)"

echo "== (5) CND $(date -Is)"
for g in $GRAPHS; do
  W=$(awk -v f="$R/cx/$g.cx" '$1==f {print $2}' $R/tods1_new3_omega.log)
  if [ -z "$W" ]; then echo "CND $g skipped: no omega"; continue; fi
  echo "CND $g sizes 2..$W $(date -Is)"
  $T python3 src-r1index/scripts/prior_sweep.py $G/$g.edges $W $g
  echo "CND $g rc=$? $(date -Is)"
done
echo "END $(date -Is)"
