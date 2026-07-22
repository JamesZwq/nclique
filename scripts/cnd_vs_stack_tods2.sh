#!/bin/bash
# THE experiment that feeds the paper (user directive 2026-07-22: every comparison must be vs CND;
# ours-vs-ours is only an ablation). Serial CND vs our full stack, same machine, same cells as
# phase-3, time AND memory, 3 trials, medians.
#
# CND rules, all learned the hard way:
#  - OMP_NUM_THREADS=1 ALWAYS. Never benchmark parallel CND (§206/§207: it is harmful on these
#    inputs and comparing to it produced two retracted claims).
#  - Run FROM THE REPO ROOT: CND opens src/nCr.txt relative to cwd; anywhere else it exits 1 on
#    every cell with "file could not be opened".
#  - 300GB address-space cap so an OOM lands as a neutral marker, not a dead box.
# Ours: the §210 full stack (SCT_DECONV + SCT_SPARSE_FP + SCT_CLAMP_PF), which phase-3 measured
# BOTH-WIN on 9/10 cells vs our own dense baseline. Its absolute numbers here should reproduce
# phase-3 within noise; the NEW information is the CND columns.
OUT=/data/wenqianz/cnd_vs_stack.out
exec >"$OUT" 2>&1
set -u
echo "=== CND VS STACK START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; WANT_HEAD="${1:-}"
cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then
  echo "FATAL: HEAD mismatch (want $WANT_HEAD)"; exit 1; fi

# --- preconditions, checked loudly ---
[ -x build/bin/degeneracy_cliques ] || { echo "FATAL: no CND binary at build/bin/degeneracy_cliques"; exit 1; }
[ -f src/nCr.txt ] || { echo "FATAL: src/nCr.txt missing (CND opens it via cwd)"; exit 1; }
echo "CND binary: $(ls -la build/bin/degeneracy_cliques | awk '{print $5, $6, $7, $8}')"
( cd region_native && g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition \
    -o /tmp/rn_stack region_native_sct_peel.cpp ); BE=$?
echo "ours build_exit=$BE"; [ $BE -eq 0 ] || exit 1

export OMP_NUM_THREADS=1
export SCT_MAX_INC=200000000
STACK="SCT_DECONV=1 SCT_SPARSE_FP=1 SCT_CLAMP_PF=1"
TMO=7200
CAP=$((300*1024*1024))   # 300GB in KB for ulimit -v
med3(){ printf '%s\n' "$1" "$2" "$3" | sort -n | sed -n 2p; }

printf '%-20s %-5s | %10s %10s | %10s %10s | %9s %9s\n' \
  GRAPH CELL CND_S CND_RSS_MB OURS_S OURS_RSS_MB TIME_RATIO MEM_RATIO
echo "----------------------------------------------------------------------------------------------"

for spec in "ca-GrQc 3 5" "ca-CondMat 3 5" "email-Eu-core 3 5" "ca-HepPh 3 5" "ca-AstroPh 4 6" \
            "com-dblp 4 6" "com-amazon.ungraph 3 5" "web-Google 3 5" "com-youtube 3 5" "web-it-2004 3 4"; do
  set -- $spec; G=$1; R=$2; S=$3
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }

  # --- CND, serial, from repo root, 3 trials ---
  c1=""; c2=""; c3=""; cfail=0
  for i in 1 2 3; do
    ( ulimit -v $CAP; cd "$REPO" && timeout $TMO /usr/bin/time -v \
        ./build/bin/degeneracy_cliques "$GR" "$R" "$S" ) >/tmp/cnd_${G}.out 2>/tmp/cnd_${G}.err
    rc=$?
    if [ $rc -ne 0 ]; then cfail=$rc; break; fi
    eval "c$i=\$(grep -oP 'Elapsed \(wall clock\).*: \K.*' /tmp/cnd_${G}.err | tail -1)"
    # normalize h:mm:ss / m:ss.xx to seconds
    eval "c$i=\$(awk -v t=\"\$c$i\" 'BEGIN{n=split(t,a,\":\"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; printf \"%.2f\", s}')"
  done
  CRSS=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/cnd_${G}.err 2>/dev/null || echo 0)

  # --- ours, full stack, 3 trials ---
  o1=""; o2=""; o3=""; ofail=0
  for i in 1 2 3; do
    timeout $TMO /usr/bin/time -v env $STACK /tmp/rn_stack "$GR" "$R" "$S" \
      >/tmp/our_${G}.out 2>/tmp/our_${G}.err || { ofail=$?; break; }
    eval "o$i=\$(grep -oP 'total=\K[0-9.]+' /tmp/our_${G}.out | tail -1)"
  done
  ORSS=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/our_${G}.err 2>/dev/null || echo 0)

  if [ $cfail -ne 0 ] && [ $ofail -ne 0 ]; then
    printf '%-20s %-5s | both failed cnd=%s ours=%s\n' "$G" "$R,$S" "$cfail" "$ofail"; continue
  elif [ $cfail -ne 0 ]; then
    om=$(med3 "$o1" "$o2" "$o3")
    printf '%-20s %-5s | %10s %10s | %10s %10d | CND-INFEASIBLE(rc=%s: 124=time,137/134=mem)\n' \
      "$G" "$R,$S" "FAIL" "-" "$om" "$((ORSS/1024))" "$cfail"; continue
  elif [ $ofail -ne 0 ]; then
    cm=$(med3 "$c1" "$c2" "$c3")
    printf '%-20s %-5s | %10s %10d | %10s %10s | OURS-FAILED(rc=%s)\n' \
      "$G" "$R,$S" "$cm" "$((CRSS/1024))" "FAIL" "-" "$ofail"; continue
  fi
  cm=$(med3 "$c1" "$c2" "$c3"); om=$(med3 "$o1" "$o2" "$o3")
  tr=$(awk -v a="$cm" -v b="$om" 'BEGIN{if(b>0)printf "%.3f",a/b; else print "-"}')
  mr=$(awk -v a="$CRSS" -v b="$ORSS" 'BEGIN{if(b>0)printf "%.3f",a/b; else print "-"}')
  printf '%-20s %-5s | %10s %10d | %10s %10d | %9s %9s\n' \
    "$G" "$R,$S" "$cm" "$((CRSS/1024))" "$om" "$((ORSS/1024))" "$tr" "$mr"
done
echo ""
echo "TIME_RATIO / MEM_RATIO = CND / ours. >1 means we win."
echo "=== CND VS STACK DONE $(date) ==="
