#!/bin/bash
# SCOUT GRID (fast-experiment protocol, 2026-07-23). Replaces the slow cnd_grid_tods2.sh.
#   ours: ONE sweep per r (SCT_SWEEP) covers all s cells -> 3 runs/graph, [nsi-cell] marginal times.
#   CND : ONE trial per cell, 1800s budget, 300GB cap; on timeout mark and SKIP larger s for that r.
#   Cells with ratio in [1/3, 3] are listed at the end as the REFINE queue (3-trial medians later).
# Usage: grid_scout_tods2.sh <expected-head> graph1 graph2 ...
OUT=/data/wenqianz/grid_scout.out
exec >"$OUT" 2>&1
set -u
echo "=== GRID SCOUT START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
WANT_HEAD="${1:-}"; shift || true
GRAPHS="${*:-web-Google com-youtube com-dblp com-amazon.ungraph web-it-2004}"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then
  echo "FATAL: HEAD mismatch"; exit 1; fi
[ -x build/bin/degeneracy_cliques ] || { echo "FATAL: no CND binary"; exit 1; }
[ -f src/nCr.txt ] || { echo "FATAL: src/nCr.txt"; exit 1; }
( cd region_native && g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition \
    -o /tmp/rn_scout region_native_sct_peel.cpp ); BE=$?
echo "ours build_exit=$BE"; [ $BE -eq 0 ] || exit 1

export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000
STACK="SCT_DECONV=1 SCT_SPARSE_FP=1 SCT_CLAMP_PF=1"
CND_TMO=1800
CAP=$((300*1024*1024))
REFINE=/data/wenqianz/grid_refine_queue.txt; : >"$REFINE"

printf '%-20s %-5s | %10s %9s | %11s %11s %9s | %9s | %s\n' \
  GRAPH CELL CND_S CND_MB OURS_MARG_S OURS_ROW_S OURS_MB T_RATIO NOTE
echo "----------------------------------------------------------------------------------------------------------"

for G in $GRAPHS; do
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  for R in 3 4 5; do
    SMAX=$((R+3))
    # ---- ours: ONE sweep covers s = R+1 .. R+3 ----
    timeout 21600 /usr/bin/time -v env $STACK SCT_SWEEP=$SMAX /tmp/rn_scout "$GR" "$R" $((R+1)) --mce-budget 7200 \
        >/tmp/s_${G}_${R}.out 2>/tmp/s_${G}_${R}.err
    ORC=$?
    ORSS=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/s_${G}_${R}.err 2>/dev/null || echo 0)
    OROW=$(grep -oP 'total=\K[0-9.]+' /tmp/s_${G}_${R}.out | tail -1)
    skipRestCND=0
    for DS in 1 2 3; do
      S=$((R+DS))
      OMARG=$(grep -oP "\[nsi-cell\] r=$R s=$S cell-time=\K[0-9.]+" /tmp/s_${G}_${R}.out | tail -1)
      [ -z "$OMARG" ] && OMARG="-"
      # ---- CND: one trial, budget, skip-on-blowup ----
      if [ $skipRestCND -eq 1 ]; then
        printf '%-20s %-5s | %10s %9s | %11s %11s %9d | %9s | CND-SKIP (smaller s already >%ss)\n' \
          "$G" "$R,$S" ">$CND_TMO" "-" "$OMARG" "${OROW:--}" "$((ORSS/1024))" "inf" "$CND_TMO"
        continue
      fi
      ( ulimit -v $CAP; cd "$REPO" && timeout $CND_TMO /usr/bin/time -v \
          ./build/bin/degeneracy_cliques "$GR" "$R" "$S" ) >/tmp/c_s.out 2>/tmp/c_s.err
      CRC=$?
      if [ $CRC -ne 0 ]; then
        skipRestCND=1
        printf '%-20s %-5s | %10s %9s | %11s %11s %9d | %9s | CND rc=%s (124=budget,others=mem)\n' \
          "$G" "$R,$S" "FAIL" "-" "$OMARG" "${OROW:--}" "$((ORSS/1024))" "inf" "$CRC"
        continue
      fi
      t=$(grep -oP 'Elapsed \(wall clock\).*: \K.*' /tmp/c_s.err | tail -1)
      CS=$(awk -v t="$t" 'BEGIN{n=split(t,a,":");s=0;for(i=1;i<=n;i++)s=s*60+a[i];printf "%.2f",s}')
      CMB=$(( $(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/c_s.err) / 1024 ))
      TR=$(awk -v a="$CS" -v b="$OMARG" 'BEGIN{if(b+0>0)printf "%.3f",a/b; else print "-"}')
      printf '%-20s %-5s | %10s %9d | %11s %11s %9d | %9s |\n' \
        "$G" "$R,$S" "$CS" "$CMB" "$OMARG" "${OROW:--}" "$((ORSS/1024))" "$TR"
      # refine queue: decision zone
      awk -v tr="$TR" -v g="$G" -v r="$R" -v s="$S" \
        'BEGIN{ if (tr+0 > 0.333 && tr+0 < 3.0) print g, r, s }' >>"$REFINE"
    done
  done
  echo ""
done
echo "=== REFINE QUEUE (ratio in [1/3,3] -- needs 3-trial medians + standalone ours) ==="
cat "$REFINE"
echo "=== GRID SCOUT DONE $(date) ==="
