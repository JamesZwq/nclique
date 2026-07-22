#!/bin/bash
# §210f: settle the sparse TIME sign on tods2 (the paper machine).
# Local (ARM) said sparse SPEEDS UP ca-AstroPh peel 1.15-1.22x; tods2 phase-2 (Xeon) said it REGRESSED
# (~0.94x). Same graph, opposite sign. This A/B compares, on tods2, peel-attributed, 3 trials:
#   dense  vs  OLD sparse (O(r^2) spGet, at commit 04a6e5f)  vs  NEW sparse (O(r) spEqVec, HEAD).
# We build BOTH binaries from git so the comparison is exact, and we read PEEL not total (the AstroPh
# MCE is tiny here, but the rule stands).
OUT=/data/wenqianz/sparse_ab_tods2.out
exec >"$OUT" 2>&1
set -u
echo "=== SPARSE A/B START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
OLD_HASH="${1:-04a6e5f}"; NEW_HASH="${2:-}"

cd "$REPO" || exit 1
git fetch -q origin
NEW_NOW=$(git rev-parse --short origin/main)
[ -n "$NEW_HASH" ] && [ "$NEW_NOW" != "$NEW_HASH" ] && { echo "note: origin/main is $NEW_NOW, wanted $NEW_HASH"; }

build_at() { # $1=hash $2=outbin
  git reset -q --hard "$1" || return 1
  echo "  build $1: $(git log --oneline -1)"
  ( cd region_native && g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition \
      -o "$2" region_native_sct_peel.cpp ) ; return $?
}
build_at "$OLD_HASH" /tmp/rn_old || { echo "FATAL old build"; exit 1; }
build_at origin/main /tmp/rn_new || { echo "FATAL new build"; exit 1; }
git reset -q --hard origin/main
grep -q 'spEqVec(sp.first, sp.second, Yscr)' region_native/region_native_sct_peel.cpp \
  && echo "guard: NEW has the O(r) spEqVec confirm" || echo "WARN: new binary may not have the fix"

export SCT_MAX_INC=200000000 OMP_NUM_THREADS=1
med3(){ printf '%s\n' "$1" "$2" "$3" | sort -n | sed -n 2p; }
peelmed(){ local BIN=$1 G=$2 R=$3 S=$4; shift 4; local p1 p2 p3
  for i in 1 2 3; do timeout 3000 env "$@" $BIN "$GD/$G.edges" "$R" "$S" >/tmp/ab.out 2>/dev/null
    eval "p$i=\$(grep -oP 'peel=\K[0-9.]+' /tmp/ab.out | tail -1)"; done
  med3 "$p1" "$p2" "$p3"; }

printf '%-14s %-5s %10s %14s %14s\n' GRAPH CELL DENSE_PEEL SPARSE_OLD SPARSE_NEW
echo "-----------------------------------------------------------------"
for spec in "ca-AstroPh 4 6" "ca-HepPh 3 5" "ca-CondMat 4 6" "com-dblp 4 6" "web-Google 3 5"; do
  set -- $spec; G=$1; R=$2; S=$3
  d=$(peelmed /tmp/rn_new "$G" "$R" "$S")                       # dense (flag off, same binary)
  o=$(peelmed /tmp/rn_old "$G" "$R" "$S" SCT_SPARSE_FP=1)       # old O(r^2)
  n=$(peelmed /tmp/rn_new "$G" "$R" "$S" SCT_SPARSE_FP=1)       # new O(r)
  so=$(awk -v a="$d" -v b="$o" 'BEGIN{if(b>0)printf "%.3f",a/b;else print"-"}')
  sn=$(awk -v a="$d" -v b="$n" 'BEGIN{if(b>0)printf "%.3f",a/b;else print"-"}')
  printf '%-14s %-5s %10s %7s(%s) %7s(%s)\n' "$G" "$R,$S" "$d" "$o" "$so" "$n" "$sn"
done
echo "=== SPARSE A/B DONE $(date) ==="
