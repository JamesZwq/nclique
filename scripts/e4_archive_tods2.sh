#!/bin/bash
# §225 E4: the QUERY axis on the acceptance roster.
#
# Compares the slim plane index (NSI3) against the MATERIALIZED ARCHIVE -- the only structure that
# answers the same queries without an index: for every r-clique, kappa in every cell of its row,
# sorted by the clique key, probed by binary search.
#
# Every methodological choice is made in the ARCHIVE's favour so the reported gap is a LOWER bound:
#   * archive rows = pattern-side r-cliques only (a strict lower bound on #r-cliques)
#   * 4-byte vertex ids, 4-byte kappa, ONE key shared across the row's cells
#   * the workload is uniform over r-cliques and hits 100% (a miss is a cheap probe)
#   * index and archive are timed by the SAME binary, on the SAME workload, back to back
#
# Timing discipline (DO_NOT_REPEAT T9): one measurement at a time, nothing else on the box.
OUT=/data/wenqianz/e4_archive.out
exec >"$OUT" 2>&1
set -u
echo "=== E4 ARCHIVE vs INDEX START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; IX=/data/wenqianz/nsi3idx
WANT_HEAD="${1:-}"; shift || true
GRAPHS="${*:-nasasrb pkustk11 pkustk13 pwtk com-amazon.ungraph com-dblp web-Google web-it-2004}"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then echo "FATAL: HEAD mismatch"; exit 1; fi
for sym in "archive-bench" "by-clique" keyBuf; do
  grep -q -- "$sym" region_native/nsi_query.cpp || { echo "FATAL: $sym missing from nsi_query.cpp"; exit 1; }
done
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -o /tmp/nq4 nsi_query.cpp || { echo "FATAL: build failed"; exit 1; }
echo "build ok"

export OMP_NUM_THREADS=1
RMIN=3; RMAX=5
NQ=20000            # per r; bench needs >=1000 and more samples tighten the median
CAPGB=200           # tods2 has 503 GB; 200 leaves ample headroom for the sort's second copy
WORK=/tmp/e4work; mkdir -p "$WORK"

for G in $GRAPHS; do
  F2="$IX/$G.nsi2"; F3="$IX/$G.nsi3"; GR="$GD/$G.edges"
  echo ""
  echo "################ $G ################"
  for f in "$F2" "$F3" "$GR"; do [ -s "$f" ] || { echo "MISSING $f"; continue 2; }; done

  echo "---- size accounting (archive vs full vs slim) ----"
  /tmp/nq4 "$F2" archive --vs "$F3"

  for R in $(seq $RMIN $RMAX); do
    Q="$WORK/$G.r$R.txt"
    /tmp/nq4 "$F2" sample $R $NQ --by-clique 2>/dev/null | awk -v r=$R 'NF{print r, $0}' > "$Q"
    n=$(wc -l < "$Q")
    if [ "$n" -lt 1000 ]; then echo "---- r=$R SKIPPED: only $n sampled cliques ----"; rm -f "$Q"; continue; fi
    echo "---- r=$R  workload=$n clique-uniform samples ----"
    echo "[INDEX  NSI3]"
    /tmp/nq4 "$F3" bench "$GR" "$Q" 2>&1 | grep -E "^benchmark|^operation|^point|^row"
    echo "[ARCHIVE]"
    /tmp/nq4 "$F2" archive-bench $R "$Q" --cap-gb $CAPGB 2>&1 | \
        grep -E "ARCHIVE-BENCH|^operation|^point|^row|MISS"
    rm -f "$Q"
  done
done
echo ""
echo "=== E4 ARCHIVE vs INDEX DONE $(date) ==="
