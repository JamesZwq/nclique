#!/bin/bash
# §236 / E7: what do the two certificates actually save?
#
# The plane build carries two transfer rules. Turn each OFF and measure what comes back:
#   SCT_NO_DIAG=1      kills the VERTICAL transfer (T5). Every row's boundary cell becomes full
#                      residue, i.e. r=4 and r=5 must peel from scratch instead of inheriting.
#   SCT_SWEEP_NOCERT=1 kills the HORIZONTAL transfer (T3). Every cell above a boundary re-derives
#                      certification per cell instead of riding the absorbing property.
#
# Both are ABLATIONS, not bugs: the engine is expected to produce the SAME kappa distribution either
# way, only slower. So the gate is the core-number distribution being identical to the default run --
# if an ablation changes an answer, the ablation is measuring something other than what it claims.
#
# Timing discipline (DO_NOT_REPEAT T9): one run at a time, machine idle.
OUT=/data/wenqianz/e7_cert.out
exec >"$OUT" 2>&1
set -u
echo "=== E7 CERTIFICATE ABLATION START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; WK=/data/wenqianz/e7; mkdir -p "$WK"
GRAPHS="${*:-nasasrb pkustk13}"
BUDGET=5400          # per run; exceeding it is itself a result

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
echo "HEAD: $(git log --oneline -1)"
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn_e7 \
    region_native_sct_peel.cpp || { echo "FATAL: build"; exit 1; }
echo "build ok"
export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000

run () {                       # run <graph> <tag> <env...>
  local G="$1" TAG="$2"; shift 2
  local f="$WK/$G.$TAG"
  /usr/bin/time -v timeout $BUDGET env "$@" SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
      /tmp/rn_e7 "$GD/$G.edges" 3 4 >"$f.out" 2>"$f.err"
  local rc=$?
  local t m
  t=$(grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K.*' "$f.err" | tail -1 | \
      awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }')
  m=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$f.err" | tail -1)
  if [ $rc -eq 124 ]; then t="OVER-${BUDGET}s"; fi
  grep -E '^core=|^\[sct-peel\]' "$f.out" | md5sum | cut -c1-12 > "$f.md5"
  printf '  %-10s rc=%-3s time=%-12s peakKB=%-12s dist=%s\n' "$TAG" "$rc" "${t:-?}" "${m:-?}" "$(cat "$f.md5")"
}

for G in $GRAPHS; do
  [ -f "$GD/$G.edges" ] || { echo "$G MISSING"; continue; }
  echo ""
  echo "################ $G ################"
  run "$G" default   DUMMY=0
  run "$G" nodiag    SCT_NO_DIAG=1
  run "$G" nocert    SCT_SWEEP_NOCERT=1
  a=$(cat "$WK/$G.default.md5" 2>/dev/null)
  for t in nodiag nocert; do
    b=$(cat "$WK/$G.$t.md5" 2>/dev/null)
    if [ -z "$b" ]; then echo "  GATE $t: NO OUTPUT (likely timed out -- that IS the result)"
    elif [ "$a" = "$b" ]; then echo "  GATE $t: distribution IDENTICAL to default"
    else echo "  GATE $t: **DISTRIBUTION DIFFERS** ($a vs $b) -- ablation is not answer-preserving"; fi
  done
  echo "---- per-cell replay times, default vs ablations ----"
  for t in default nodiag nocert; do
    echo "  [$t]"
    grep -E 'nsi-plane-cell' "$WK/$G.$t.out" 2>/dev/null | \
      sed 's/\[nsi-plane-cell\] /    /' | cut -c1-110
  done
done
echo ""
echo "=== E7 DONE $(date) ==="
