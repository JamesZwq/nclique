#!/bin/bash
# §228 A/B on the acceptance roster: the plane index build, before vs after moving the diagonal
# certificate ahead of the leaf maps (and dropping the per-pattern host list / duplicate comp key).
#
# The gate is the written index being BYTE-IDENTICAL, which is the strongest available: if a single
# certification decision or cP changed, the file changes. Time and peak RSS are both reported because
# both must improve -- a build that trades memory for speed does not count here.
#
# Timing discipline (DO_NOT_REPEAT T9): one build at a time, nothing else on the box.
OUT=/data/wenqianz/plane228.out
exec >"$OUT" 2>&1
set -u
echo "=== PLANE 228 A/B START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; WK=/data/wenqianz/p228; mkdir -p "$WK"
WANT_HEAD="${1:-}"; BASE="${2:-e176118}"; shift 2 2>/dev/null || true
GRAPHS="${*:-nasasrb pkustk11 pwtk com-amazon.ungraph com-dblp web-it-2004 pkustk13 web-Google}"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"; echo "BASE: $(git log --oneline -1 $BASE)"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then echo "FATAL: HEAD mismatch"; exit 1; fi
grep -q "needLeafMaps" region_native/region_native_sct_peel.cpp || { echo "FATAL: 228 not in tree"; exit 1; }
git show "$BASE:region_native/region_native_sct_peel.cpp" > /tmp/rn_base.cpp || exit 1
grep -q "needLeafMaps" /tmp/rn_base.cpp && { echo "FATAL: BASE already has 228"; exit 1; }
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn228_base /tmp/rn_base.cpp || { echo "FATAL: base build"; exit 1; }
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn228_opt region_native_sct_peel.cpp || { echo "FATAL: opt build"; exit 1; }
echo "builds ok"

export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000
RMIN=3; RMAX=5; SMAX=8

printf '%-22s %10s %10s %8s %12s %12s %8s %s\n' GRAPH BASE_S OPT_S SPEEDUP BASE_KB OPT_KB MEMGAIN GATE
echo "----------------------------------------------------------------------------------------------------"
for G in $GRAPHS; do
  GR="$GD/$G.edges"; [ -f "$GR" ] || { printf '%-22s MISSING\n' "$G"; continue; }
  R=()
  for V in base opt; do
    F="$WK/$G.$V.nsi2"; rm -f "$F"
    /usr/bin/time -v env SCT_RSWEEP=1 SCT_RMIN=$RMIN SCT_RMAX=$RMAX SCT_SMAX=$SMAX \
        SCT_INDEX_OUT="$F" /tmp/rn228_$V "$GR" $RMIN $((RMIN+1)) \
        >"$WK/$G.$V.out" 2>"$WK/$G.$V.err"
    rc=$?
    if [ $rc -ne 0 ] || [ ! -s "$F" ]; then printf '%-22s %s BUILD FAILED rc=%s\n' "$G" "$V" "$rc"; R=(); break; fi
    t=$(grep -oP 'Elapsed \(wall clock\) time.*?: \K.*' "$WK/$G.$V.err" | tail -1 | \
        awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }')
    m=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$WK/$G.$V.err" | tail -1)
    R+=("$t" "$m")
  done
  [ ${#R[@]} -eq 4 ] || continue
  if cmp -s "$WK/$G.base.nsi2" "$WK/$G.opt.nsi2"; then GATE="BYTE-IDENTICAL"; else GATE="**INDEX DIFFERS**"; fi
  # DO_NOT_REPEAT T11: the server awk rejects a ternary inside a printf ARGUMENT.
  SP=$(awk -v a="${R[0]}" -v b="${R[2]}" 'BEGIN{ if (b>0) printf "%.2f", a/b; else printf "0" }')
  MG=$(awk -v a="${R[1]}" -v b="${R[3]}" 'BEGIN{ if (b>0) printf "%.2f", a/b; else printf "0" }')
  printf '%-22s %10s %10s %7sx %12s %12s %7sx %s\n' "$G" "${R[0]}" "${R[2]}" "$SP" "${R[1]}" "${R[3]}" "$MG" "$GATE"
  rm -f "$WK/$G.base.nsi2" "$WK/$G.opt.nsi2"
done

echo ""
echo "=== per-column phase split (opt runs): where the time went, and what was skipped ==="
for G in $GRAPHS; do
  [ -f "$WK/$G.opt.out" ] || continue
  echo "---- $G ----"
  grep -E 'leaf-maps=|phase-split|complete r=' "$WK/$G.opt.out"
done
echo ""
echo "=== PLANE 228 A/B DONE $(date) ==="
