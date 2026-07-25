#!/bin/bash
# §223 NSI3 on the acceptance roster: build both the full (NSI2) and slim (NSI3) plane indexes,
# gate the slim one answer-for-answer against the full one, then report size / load / query.
# Gating is per r with r-cliques SAMPLED FOR THAT r: feeding one r's cliques to another column
# tests non-cliques and produces meaningless verdicts (§223 lesson).
OUT=/data/wenqianz/nsi3.out
exec >"$OUT" 2>&1
set -u
echo "=== NSI3 ROSTER START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; IX=/data/wenqianz/nsi3idx
WANT_HEAD="${1:-}"; shift || true
GRAPHS="${*:-nasasrb pkustk11 pkustk13 pwtk com-amazon.ungraph com-dblp web-Google web-it-2004}"
mkdir -p "$IX"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then echo "FATAL: HEAD mismatch"; exit 1; fi
for sym in SCT_INDEX_SLIM NSI3 cpFromProfiles; do
  grep -q "$sym" region_native/region_native_sct_peel.cpp region_native/nsi_query.cpp || { echo "FATAL: $sym missing"; exit 1; }
done
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn3 region_native_sct_peel.cpp || exit 1
g++ -O3 -march=native -std=c++17 -o /tmp/nq3 nsi_query.cpp || exit 1
echo "builds ok"

export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000
RMIN=3; RMAX=5; SMAX=8
NQ=50000

printf '%-22s %14s %14s %8s %12s %12s %10s %s\n' GRAPH NSI2_B NSI3_B SHRINK NSI2_B/RC NSI3_B/RC LOAD_MS GATE
echo "--------------------------------------------------------------------------------------------------------"

for G in $GRAPHS; do
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  F2="$IX/$G.nsi2"; F3="$IX/$G.nsi3"
  timeout 21600 env SCT_RSWEEP=1 SCT_RMIN=$RMIN SCT_RMAX=$RMAX SCT_SMAX=$SMAX SCT_INDEX_OUT="$F2" \
      /tmp/rn3 "$GR" $RMIN $((RMIN+1)) >/tmp/n2_$G.out 2>/tmp/n2_$G.err
  R2=$?
  timeout 21600 env SCT_RSWEEP=1 SCT_RMIN=$RMIN SCT_RMAX=$RMAX SCT_SMAX=$SMAX SCT_INDEX_OUT="$F3" \
      SCT_INDEX_SLIM=1 /tmp/rn3 "$GR" $RMIN $((RMIN+1)) >/tmp/n3_$G.out 2>/tmp/n3_$G.err
  R3=$?
  if [ $R2 -ne 0 ] || [ $R3 -ne 0 ] || [ ! -s "$F2" ] || [ ! -s "$F3" ]; then
    printf '%-22s BUILD FAILED (nsi2 rc=%s nsi3 rc=%s)\n' "$G" "$R2" "$R3"; continue
  fi
  B2=$(stat -c%s "$F2"); B3=$(stat -c%s "$F3")

  GATE=PASS; ROWS=0
  for R in $(seq $RMIN $RMAX); do
    /tmp/nq3 "$F2" sample $R $NQ >/tmp/q_$G$R.txt 2>/dev/null
    n=$(wc -l </tmp/q_$G$R.txt)
    [ "$n" -lt 100 ] && { GATE="NO-SAMPLE(r=$R)"; continue; }
    /tmp/nq3 "$F2" rowfile $R /tmp/q_$G$R.txt >/tmp/a2.txt 2>/dev/null
    /tmp/nq3 "$F3" rowfile $R /tmp/q_$G$R.txt >/tmp/a3.txt 2>/dev/null
    if [ ! -s /tmp/a2.txt ]; then GATE="NO-REF(r=$R)"; continue; fi
    if cmp -s /tmp/a2.txt /tmp/a3.txt; then ROWS=$((ROWS+n)); else GATE="**DIFFERS(r=$R)**"; fi
  done

  L3=$(/tmp/nq3 "$F3" stats 2>/dev/null | grep -oP 'load=\K[0-9.]+' | head -1)
  RC=$(grep -oP 'r-cliques=\K[0-9]+' /tmp/n2_$G.out | awk '{t+=$1} END{print t+0}')
  [ "${RC:-0}" -eq 0 ] && RC=1
  printf '%-22s %14s %14s %7sx %12s %12s %10s %s(%s rows)\n' \
    "$G" "$B2" "$B3" "$(awk -v a=$B2 -v b=$B3 'BEGIN{printf "%.1f",a/b}')" \
    "$(awk -v a=$B2 -v r=$RC 'BEGIN{printf "%.4f",a/r}')" \
    "$(awk -v a=$B3 -v r=$RC 'BEGIN{printf "%.4f",a/r}')" "${L3:--}" "$GATE" "$ROWS"

  { echo "### $G"; grep -E 'nsi3-slim|nsi2-bytes' /tmp/n3_$G.err /tmp/n2_$G.err 2>/dev/null | sed 's/^/    /'; } \
    >>/data/wenqianz/nsi3_detail.out
done
echo ""
echo "=== NSI3 ROSTER DONE $(date) ==="
