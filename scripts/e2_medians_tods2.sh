#!/bin/bash
# E2: 3-trial medians for every number the paper quotes.
#
# Design rules (DO_NOT_REPEAT T9/T10/T12 + benchmark-setup rules):
#   * strictly one measurement at a time, machine to ourselves
#   * resume-friendly: results land in a CSV; rows already OK are skipped on restart
#   * failures are censored observations: established single-run failures (the 30
#     infeasible grid cells) are NOT re-run x3; the one >5400s ablation gets 1 trial
#   * CND runs from the repo root (src/nCr.txt), serial (OMP_NUM_THREADS=1)
#   * every run wrapped in /usr/bin/time -v; full stdout+stderr kept per trial
#
# Order: A plane-build x3 (8 graphs) -> B ablations x3 (NS, PK3) ->
#        C grid ours-sweep x3 + CND feasible cells x3 -> D WI(3,4) standalone x3 ->
#        E query bench x3.  Estimated total ~17 h serial.
OUT=/data/wenqianz/e2
CSV=$OUT/results.csv
LOG=$OUT/e2_driver.out
mkdir -p $OUT/logs
exec >>"$LOG" 2>&1
set -u
echo "=== E2 MEDIANS START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
echo "HEAD: $(git log --oneline -1)"
grep -q "needLeafMaps" region_native/region_native_sct_peel.cpp || { echo "FATAL: engine symbol missing"; exit 1; }
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn_e2 region_native_sct_peel.cpp || { echo "FATAL build rn"; exit 1; }
g++ -O3 -march=native -std=c++17 -o /tmp/nq_e2 nsi_query.cpp || { echo "FATAL build nq"; exit 1; }
cd "$REPO"
[ -x build/bin/degeneracy_cliques ] || { echo "FATAL: CND binary missing"; exit 1; }
[ -f src/nCr.txt ] || cp pivoter/nCr.txt src/nCr.txt || { echo "FATAL: nCr.txt"; exit 1; }
[ -f "$CSV" ] || echo "exp,graph,config,trial,status,secs,peakKB" > "$CSV"
export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000

have () { grep -q "^$1,$2,$3,$4,OK" "$CSV"; }
tparse () { grep -oP 'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): \K.*' "$1" | tail -1 | \
  awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else if (NF==2) print $1*60+$2; else print $1 }'; }
mparse () { grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' "$1" | tail -1; }

run1 () {  # run1 EXP GRAPH CONFIG TRIAL BUDGET CMD...
  local E=$1 G=$2 C=$3 T=$4 B=$5; shift 5
  have "$E" "$G" "$C" "$T" && { echo "skip $E,$G,$C,$T"; return; }
  local f=$OUT/logs/${E}_${G}_${C}_t${T}
  /usr/bin/time -v timeout "$B" "$@" >"$f.out" 2>"$f.err"
  local rc=$?
  local secs mem st
  secs=$(tparse "$f.err"); mem=$(mparse "$f.err")
  if [ $rc -eq 0 ]; then st=OK; elif [ $rc -eq 124 ]; then st=TIMEOUT; else st=ERR$rc; fi
  echo "$E,$G,$C,$T,$st,${secs:-},${mem:-}" >> "$CSV"
  echo "$(date +%H:%M) $E,$G,$C,$T -> $st ${secs:-}s"
}

GRAPHS="web-it-2004 web-Google com-dblp com-amazon.ungraph nasasrb pkustk11 pkustk13 pwtk"

# ---------- A. whole-plane build x3 ----------
for T in 1 2 3; do
  for G in $GRAPHS; do
    run1 plane "$G" default "$T" 7200 env SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
        /tmp/rn_e2 "$GD/$G.edges" 3 4
  done
done

# ---------- B. ablations x3 (censored config: 1 trial) ----------
for T in 1 2 3; do
  for G in nasasrb pkustk13; do
    run1 abl "$G" nodiag "$T" 7200 env SCT_NO_DIAG=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
        /tmp/rn_e2 "$GD/$G.edges" 3 4
  done
  run1 abl nasasrb nocert "$T" 7200 env SCT_SWEEP_NOCERT=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
      /tmp/rn_e2 "$GD/nasasrb.edges" 3 4
done
run1 abl pkustk13 nocert 1 5400 env SCT_SWEEP_NOCERT=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
    /tmp/rn_e2 "$GD/pkustk13.edges" 3 4

# ---------- C. grid: ours sweep x3 per (graph, r); CND x3 per FEASIBLE cell ----------
for T in 1 2 3; do
  for G in $GRAPHS; do
    for R in 3 4 5; do
      run1 sweep "$G" "r$R" "$T" 7200 env SCT_SWEEP=1 SCT_SMAX=$((R+3>8?8:R+3)) \
          /tmp/rn_e2 "$GD/$G.edges" $R $((R+1))
    done
  done
done
# feasible CND cells only (established failures are censored; not re-run)
CND_CELLS="web-Google:3,4 web-Google:3,5 web-Google:3,6 web-Google:4,5 web-Google:4,6 web-Google:4,7 web-Google:5,6 web-Google:5,7 web-Google:5,8 \
com-dblp:3,4 com-dblp:3,5 com-dblp:3,6 com-dblp:4,5 com-dblp:4,6 com-dblp:4,7 com-dblp:5,6 com-dblp:5,7 com-dblp:5,8 \
com-amazon.ungraph:3,4 com-amazon.ungraph:3,5 com-amazon.ungraph:3,6 com-amazon.ungraph:4,5 com-amazon.ungraph:4,6 com-amazon.ungraph:4,7 com-amazon.ungraph:5,6 com-amazon.ungraph:5,7 com-amazon.ungraph:5,8 \
nasasrb:3,4 nasasrb:3,5 nasasrb:3,6 nasasrb:4,5 nasasrb:4,6 nasasrb:4,7 \
pkustk11:3,4 pkustk11:3,5 pkustk11:3,6 \
pkustk13:3,4 pkustk13:3,5 pkustk13:3,6 \
pwtk:3,4 pwtk:3,5 pwtk:3,6"
for T in 1 2 3; do
  for GC in $CND_CELLS; do
    G="${GC%%:*}"; RS="${GC##*:}"; R="${RS%%,*}"; S="${RS##*,}"
    run1 cnd "$G" "$R-$S" "$T" 1800 ./build/bin/degeneracy_cliques "$GD/$G.edges" "$R" "$S"
  done
done

# ---------- D. WI (3,4) standalone x3, both sides ----------
for T in 1 2 3; do
  run1 wi34 web-it-2004 ours "$T" 7200 env /tmp/rn_e2 "$GD/web-it-2004.edges" 3 4
  run1 wi34 web-it-2004 cnd "$T" 7200 ./build/bin/degeneracy_cliques "$GD/web-it-2004.edges" 3 4
done

# ---------- E. query bench x3 per graph ----------
IX=/data/wenqianz/nsi3idx
for T in 1 2 3; do
  for G in $GRAPHS; do
    [ -s "$IX/$G.nsi4" ] || continue
    for R in 3 4 5; do
      Q=/tmp/e2q_$G$R.txt
      if [ -s "$IX/$G.nsi2" ]; then /tmp/nq_e2 "$IX/$G.nsi2" sample $R 20000 --by-clique 2>/dev/null | awk -v r=$R 'NF{print r,$0}' >"$Q"; fi
      [ -s "$Q" ] || continue
      if have query "$G" "r$R-warmns" "$T"; then echo "skip query,$G,r$R,$T"; rm -f "$Q"; continue; fi
      f=$OUT/logs/query_${G}_r${R}_t${T}
      timeout 1200 /tmp/nq_e2 "$IX/$G.nsi4" bench "$GD/$G.edges" "$Q" >"$f.out" 2>"$f.err"
      rc=$?
      w=$(awk '/^point +kernel/{print $3}' "$f.out"); c=$(awk '/^point +kernel/{print $5}' "$f.out")
      echo "query,$G,r$R-warmns,$T,$([ $rc -eq 0 ] && echo OK || echo ERR$rc),${w:-},${c:-}" >> "$CSV"
      rm -f "$Q"
    done
  done
done

echo "=== E2 MEDIANS DONE $(date) ==="
