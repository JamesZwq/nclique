#!/bin/bash
# §227 NSI4 on the acceptance roster: re-encode each existing NSI3 as NSI4, gate it answer-for-answer
# against BOTH the slim NSI3 and the full NSI2, and report the size / load table.
#
# NSI4 is an ENCODING of NSI3, not a second design: it decodes into the same in-memory structures and
# runs the same query kernel. So the only thing that can go wrong is the encoding, and the only gate
# that can catch it is answer equality (DO_NOT_REPEAT T6: a size-only or structure-only check once
# called a silently-all-zero index a triumph).
#
# Converting instead of rebuilding is deliberate: the engine writes NSI4 directly, and its output is
# gated byte-for-byte against the converter locally, so a 40-minute roster rebuild buys nothing here.
OUT=/data/wenqianz/nsi4.out
exec >"$OUT" 2>&1
set -u
echo "=== NSI4 ROSTER START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; IX=/data/wenqianz/nsi3idx
WANT_HEAD="${1:-}"; shift || true
GRAPHS="${*:-nasasrb pkustk11 pkustk13 pwtk com-amazon.ungraph com-dblp web-Google web-it-2004}"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then echo "FATAL: HEAD mismatch"; exit 1; fi
for sym in loadNSI4 NSI4 "mode == \"pack\"" "mode == \"anatomy\""; do
  grep -q -- "$sym" region_native/nsi_query.cpp || { echo "FATAL: $sym missing"; exit 1; }
done
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -o /tmp/nq5 nsi_query.cpp || { echo "FATAL: build failed"; exit 1; }
echo "build ok"

NQ=20000
printf '%-22s %12s %12s %12s %8s %8s %10s %10s %s\n' \
       GRAPH NSI2_B NSI3_B NSI4_B 3to4 2to4 LOAD3_MS LOAD4_MS GATE
echo "---------------------------------------------------------------------------------------------------------------"

for G in $GRAPHS; do
  F2="$IX/$G.nsi2"; F3="$IX/$G.nsi3"; F4="$IX/$G.nsi4"
  [ -s "$F3" ] || { printf '%-22s MISSING %s\n' "$G" "$F3"; continue; }

  /tmp/nq5 "$F3" pack "$F4" > /tmp/pack_$G.txt 2>&1 || { printf '%-22s PACK FAILED\n' "$G"; cat /tmp/pack_$G.txt; continue; }
  B2=$([ -s "$F2" ] && stat -c%s "$F2" || echo 0); B3=$(stat -c%s "$F3"); B4=$(stat -c%s "$F4")

  GATE=PASS; ROWS=0
  for R in 3 4 5; do
    Q=/tmp/s_$G$R.txt
    if [ -s "$F2" ]; then /tmp/nq5 "$F2" sample $R $NQ --by-clique >"$Q" 2>/dev/null
    else                  /tmp/nq5 "$F3" sample $R $NQ            >"$Q" 2>/dev/null; fi
    n=$(wc -l <"$Q")
    if [ "$n" -lt 100 ]; then GATE="NO-SAMPLE(r=$R)"; rm -f "$Q"; continue; fi
    /tmp/nq5 "$F3" rowfile $R "$Q" >/tmp/g3.txt 2>/dev/null
    /tmp/nq5 "$F4" rowfile $R "$Q" >/tmp/g4.txt 2>/dev/null
    if [ ! -s /tmp/g3.txt ]; then GATE="NO-REF(r=$R)"; rm -f "$Q"; continue; fi
    if ! cmp -s /tmp/g3.txt /tmp/g4.txt; then GATE="**DIFFERS-vs-NSI3(r=$R)**"; rm -f "$Q"; continue; fi
    if [ -s "$F2" ]; then
      /tmp/nq5 "$F2" rowfile $R "$Q" >/tmp/g2.txt 2>/dev/null
      cmp -s /tmp/g2.txt /tmp/g4.txt || GATE="**DIFFERS-vs-NSI2(r=$R)**"
    fi
    ROWS=$((ROWS+n)); rm -f "$Q"
  done
  [ "$GATE" = PASS ] && GATE="PASS($ROWS rows)"

  L3=$(/tmp/nq5 "$F3" stats 2>/dev/null | grep -oP 'load=\K[0-9.]+' | head -1)
  L4=$(/tmp/nq5 "$F4" stats 2>/dev/null | grep -oP 'load=\K[0-9.]+' | head -1)
  R34=$(awk -v a="$B3" -v b="$B4" 'BEGIN{printf "%.2f", b>0?a/b:0}')
  R24=$(awk -v a="$B2" -v b="$B4" 'BEGIN{printf "%.1f", b>0?a/b:0}')
  printf '%-22s %12s %12s %12s %7sx %7sx %10s %10s %s\n' "$G" "$B2" "$B3" "$B4" "$R34" "$R24" "${L3:-?}" "${L4:-?}" "$GATE"
done

echo ""
echo "=== §226 byte anatomy of what NSI3 still keeps (the reason NSI4 exists) ==="
for G in $GRAPHS; do
  F3="$IX/$G.nsi3"; [ -s "$F3" ] || continue
  echo "---- $G ----"
  /tmp/nq5 "$F3" anatomy 2>/dev/null
done
echo ""
echo "=== NSI4 ROSTER DONE $(date) ==="
