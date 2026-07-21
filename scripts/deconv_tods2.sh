#!/bin/bash
# §208 G2 cross-graph validation on tods2.
# Baseline (production supFast) vs SCT_DECONV=1, SAME machine, SERIAL, one cell at a time.
# Purpose: a single cell CANNOT decide this. Different graphs put the cost in different
# places (ca-HepPh supInit-bound, ca-AstroPh addDelta-bound, web-it MCE-bound), so the
# deconvolution is expected to help a lot / a little / not at all across this set, and may
# LOSE where M is small. We need the whole spread before any conclusion.
OUT=/data/wenqianz/deconv_tods2.out
exec >"$OUT" 2>&1
set -u
echo "=== DECONV CROSS-GRAPH START $(date) ==="

REPO=/home/wenqianz/UNSW/pivoter
GD=/data/wenqianz
WANT_HEAD="${1:-}"

cd "$REPO" || { echo "FATAL: no repo"; exit 1; }
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1)
echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then
  echo "FATAL: HEAD does not contain expected hash $WANT_HEAD -- refusing to run"; exit 1
fi

# symbol guards: proves we built the source we think we did (§202 warning 1)
echo "guard supDeconv   = $(grep -c 'supDeconv' region_native/region_native_sct_peel.cpp)"
echo "guard SCT_DECONV  = $(grep -c 'SCT_DECONV' region_native/region_native_sct_peel.cpp)"
echo "guard G1 probe    = $(grep -c '208 G1' region_native/region_native_sct_peel.cpp)"
if [ "$(grep -c 'supDeconv' region_native/region_native_sct_peel.cpp)" -lt 3 ]; then
  echo "FATAL: supDeconv missing from source -- refusing to run"; exit 1
fi

cd "$REPO/region_native" || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition \
    -o /tmp/rn_dec region_native_sct_peel.cpp
BE=$?
echo "build_exit=$BE"
[ $BE -eq 0 ] || { echo "FATAL: build failed"; exit 1; }
echo "binary=$(ls -la /tmp/rn_dec | awk '{print $5}')"

export SCT_MAX_INC=200000000
export OMP_NUM_THREADS=1
TMO=2400

# cheap-first so partial results are useful if the whole thing is cut short
CELLS="
ca-GrQc 3 5
ca-GrQc 4 6
email-Eu-core 3 5
ca-CondMat 3 5
ca-CondMat 4 6
ca-HepPh 3 5
ca-AstroPh 4 6
com-amazon.ungraph 3 5
com-dblp 4 6
web-Stanford 3 4
web-it-2004 3 4
web-Google 3 5
com-youtube 3 5
"

printf '%-22s %-5s %10s %10s %8s %10s %10s %8s %s\n' \
  GRAPH CELL BASE_S DEC_S SPEEDUP BASE_RSS_MB DEC_RSS_MB MEM_RATIO CORE
echo "-----------------------------------------------------------------------------------------------"

echo "$CELLS" | while read -r G R S; do
  [ -z "${G:-}" ] && continue
  GR="$GD/$G.edges"
  if [ ! -f "$GR" ]; then echo "$G ($R,$S) GRAPH_MISSING"; continue; fi

  timeout $TMO /usr/bin/time -v /tmp/rn_dec "$GR" "$R" "$S" \
      >/tmp/b_$G$R$S.out 2>/tmp/b_$G$R$S.err
  BX=$?
  timeout $TMO /usr/bin/time -v env SCT_DECONV=1 /tmp/rn_dec "$GR" "$R" "$S" \
      >/tmp/d_$G$R$S.out 2>/tmp/d_$G$R$S.err
  DX=$?

  if [ $BX -ne 0 ] || [ $DX -ne 0 ]; then
    echo "$G ($R,$S) NONZERO_EXIT base=$BX dec=$DX  (124=timeout ${TMO}s, 7=SCT_MAX_INC guard)"
    continue
  fi

  bt=$(grep -oP 'total=\K[0-9.]+' /tmp/b_$G$R$S.out | tail -1)
  dt=$(grep -oP 'total=\K[0-9.]+' /tmp/d_$G$R$S.out | tail -1)
  br=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/b_$G$R$S.err)
  dr=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/d_$G$R$S.err)
  grep -E '^core=' /tmp/b_$G$R$S.out >/tmp/b_$G$R$S.core
  grep -E '^core=' /tmp/d_$G$R$S.out >/tmp/d_$G$R$S.core
  nb=$(wc -l </tmp/b_$G$R$S.core)

  if [ "$nb" -eq 0 ]; then CORE="NO_BASELINE_OUTPUT"
  elif diff -q /tmp/b_$G$R$S.core /tmp/d_$G$R$S.core >/dev/null; then CORE="BITEXACT($nb)"
  else CORE="**DIFFERS**"; fi

  sp=$(awk -v a="$bt" -v b="$dt" 'BEGIN{if(b>0)printf "%.3f", a/b; else print "NA"}')
  mr=$(awk -v a="$br" -v b="$dr" 'BEGIN{if(a>0)printf "%.3f", b/a; else print "NA"}')
  printf '%-22s %-5s %10s %10s %8s %10d %10d %8s %s\n' \
    "$G" "$R,$S" "$bt" "$dt" "$sp" "$((br/1024))" "$((dr/1024))" "$mr" "$CORE"

  # per-cell detail for the record
  {
    echo "### $G ($R,$S)"
    grep -E 'TIMING|supinit-prof|deconv .208|peel-prof .120. segments' /tmp/b_$G$R$S.out /tmp/b_$G$R$S.err
    grep -E 'TIMING|supinit-prof|deconv .208|peel-prof .120. segments' /tmp/d_$G$R$S.out /tmp/d_$G$R$S.err
  } >>/data/wenqianz/deconv_tods2_detail.out
done

echo ""
echo "=== DECONV CROSS-GRAPH DONE $(date) ==="
