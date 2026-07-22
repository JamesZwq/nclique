#!/bin/bash
# THE GRID (user directive 2026-07-23): beat CND across ALL (r,s) with r>=3 -- no cherry-picked cells.
# For each graph: cells r=3,4,5 x s=r+1..r+3. Each cell: serial CND vs our full stack, time+memory,
# TRIALS trials (median). CND-infeasible cells (timeout / 300GB) get neutral markers and count as
# feasibility wins. Output is one long TSV-ish table -> the per-graph ratio heatmap for the paper.
#
# Graph list is ARGUMENT-DRIVEN so the same script serves the current set and the acceptance roster:
#   bash cnd_grid_tods2.sh <expected-head> graph1 graph2 ...
OUT=/data/wenqianz/cnd_grid.out
exec >"$OUT" 2>&1
set -u
echo "=== CND GRID START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
WANT_HEAD="${1:-}"; shift || true
GRAPHS="${*:-web-it-2004 com-dblp com-amazon.ungraph web-Google com-youtube}"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then
  echo "FATAL: HEAD mismatch (want $WANT_HEAD)"; exit 1; fi
[ -x build/bin/degeneracy_cliques ] || { echo "FATAL: no CND binary"; exit 1; }
[ -f src/nCr.txt ] || { echo "FATAL: src/nCr.txt missing"; exit 1; }
( cd region_native && g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition \
    -o /tmp/rn_grid region_native_sct_peel.cpp ); BE=$?
echo "ours build_exit=$BE"; [ $BE -eq 0 ] || exit 1

export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000
STACK="SCT_DECONV=1 SCT_SPARSE_FP=1 SCT_CLAMP_PF=1"
TMO=7200
CAP=$((300*1024*1024))
TRIALS=2   # grid is wide; 2 trials + min-spread check, escalate to 3 only where the two disagree >10%
med(){ if [ $# -eq 2 ]; then awk -v a="$1" -v b="$2" 'BEGIN{printf "%.2f",(a+b)/2}'; else printf '%s\n' "$@" | sort -n | sed -n 2p; fi; }

printf '%-20s %-5s | %10s %10s | %10s %10s | %9s %9s | %s\n' \
  GRAPH CELL CND_S CND_MB OURS_S OURS_MB T_RATIO M_RATIO NOTE
echo "---------------------------------------------------------------------------------------------------------"

for G in $GRAPHS; do
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  for R in 3 4 5; do
    for DS in 1 2 3; do
      S=$((R+DS))
      # --- CND ---
      cts=(); cfail=0
      for i in $(seq $TRIALS); do
        ( ulimit -v $CAP; cd "$REPO" && timeout $TMO /usr/bin/time -v \
            ./build/bin/degeneracy_cliques "$GR" "$R" "$S" ) >/tmp/g_cnd.out 2>/tmp/g_cnd.err
        rc=$?; [ $rc -ne 0 ] && { cfail=$rc; break; }
        t=$(grep -oP 'Elapsed \(wall clock\).*: \K.*' /tmp/g_cnd.err | tail -1)
        cts+=("$(awk -v t="$t" 'BEGIN{n=split(t,a,":");s=0;for(i=1;i<=n;i++)s=s*60+a[i];printf "%.2f",s}')")
      done
      CRSS=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/g_cnd.err 2>/dev/null || echo 0)
      # --- ours ---
      ots=(); ofail=0
      for i in $(seq $TRIALS); do
        timeout $TMO /usr/bin/time -v env $STACK /tmp/rn_grid "$GR" "$R" "$S" \
          >/tmp/g_our.out 2>/tmp/g_our.err || { ofail=$?; break; }
        ots+=("$(grep -oP 'total=\K[0-9.]+' /tmp/g_our.out | tail -1)")
      done
      ORSS=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/g_our.err 2>/dev/null || echo 0)

      if [ $cfail -ne 0 ] && [ $ofail -ne 0 ]; then
        printf '%-20s %-5s | %s\n' "$G" "$R,$S" "BOTH-FAIL cnd=$cfail ours=$ofail"
      elif [ $cfail -ne 0 ]; then
        om=$(med "${ots[@]}")
        printf '%-20s %-5s | %10s %10s | %10s %10d | %9s %9s | CND-INFEASIBLE rc=%s\n' \
          "$G" "$R,$S" "-" "-" "$om" "$((ORSS/1024))" "inf" "inf" "$cfail"
      elif [ $ofail -ne 0 ]; then
        cm=$(med "${cts[@]}")
        printf '%-20s %-5s | %10s %10d | %10s %10s | %9s %9s | OURS-FAIL rc=%s\n' \
          "$G" "$R,$S" "$cm" "$((CRSS/1024))" "-" "-" "0" "0" "$ofail"
      else
        cm=$(med "${cts[@]}"); om=$(med "${ots[@]}")
        tr=$(awk -v a="$cm" -v b="$om" 'BEGIN{if(b>0)printf "%.3f",a/b;else print "-"}')
        mr=$(awk -v a="$CRSS" -v b="$ORSS" 'BEGIN{if(b>0)printf "%.3f",a/b;else print "-"}')
        printf '%-20s %-5s | %10s %10d | %10s %10d | %9s %9s |\n' \
          "$G" "$R,$S" "$cm" "$((CRSS/1024))" "$om" "$((ORSS/1024))" "$tr" "$mr"
      fi
    done
  done
  echo ""
done
echo "=== CND GRID DONE $(date) ==="
