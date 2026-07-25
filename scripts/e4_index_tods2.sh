#!/bin/bash
# E4: THE INDEX EXPERIMENT (the missing leg of the index-first paper).
# Per roster graph x r: build the NSI1 index over a whole s-row, then compare against the
# MATERIALIZED ARCHIVE (kappa per r-clique per cell, sorted, binary-searched) on both axes:
#   SIZE  index bytes  vs  archive bytes
#   QUERY nsi_query kappa()  vs  binary search over the archive, SAME sampled workload
# Where the archive exceeds the cap it CANNOT be built: that is the headline regime, reported as such.
# Index build runs on the SWEEP path, which carries the full §210 stack (deconv+sparse+clamp-pf).
OUT=/data/wenqianz/e4_index.out
exec >"$OUT" 2>&1
set -u
echo "=== E4 INDEX START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; ID=/data/wenqianz/e4idx
WANT_HEAD="${1:-}"; shift || true
GRAPHS="${*:-nasasrb pkustk11 pkustk13 pwtk com-amazon.ungraph com-dblp web-Google web-it-2004}"
mkdir -p "$ID"

cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then echo "FATAL: HEAD mismatch"; exit 1; fi
for sym in SCT_INDEX_OUT nsi_baseline.cpp; do
  [ -e "region_native/$sym" ] || grep -q "$sym" region_native/region_native_sct_peel.cpp || { echo "FATAL: $sym missing"; exit 1; }
done
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn_e4 region_native_sct_peel.cpp || exit 1
g++ -O3 -march=native -std=c++17 -o /tmp/nsi_query nsi_query.cpp || exit 1
g++ -O3 -march=native -std=c++17 -o /tmp/nsi_base nsi_baseline.cpp || exit 1
echo "builds ok"

export OMP_NUM_THREADS=1 SCT_MAX_INC=500000000
STACK="SCT_DECONV=1 SCT_SPARSE_FP=1 SCT_CLAMP_PF=1"
NQ=200000
CAP=4e8

printf '%-20s %-4s %-10s %10s %12s %10s %12s %10s %10s %10s\n' \
  GRAPH r CELLS BUILD_S IDX_MB IDX_B/RC ARCH_MB ARCH_B/RC SIZE_X OURS_NS
echo "--------------------------------------------------------------------------------------------------------------"

for G in $GRAPHS; do
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  for R in 3 4; do
    S0=$((R+1)); SMAX=$((R+3))
    IDXF="$ID/${G}_r${R}.nsi"
    T0=$(date +%s)
    timeout 14400 env $STACK SCT_SWEEP=$SMAX SCT_INDEX_OUT="$IDXF" \
        /tmp/rn_e4 "$GR" "$R" "$S0" >/tmp/e4_${G}_${R}.out 2>/tmp/e4_${G}_${R}.err
    RC=$?; T1=$(date +%s)
    if [ $RC -ne 0 ] || [ ! -s "$IDXF" ]; then
      printf '%-20s %-4s %-10s %s\n' "$G" "$R" "$S0..$SMAX" "BUILD FAILED rc=$RC ($(($T1-$T0))s)"
      continue
    fi
    BUILD=$((T1-T0))
    # sample a workload of REAL r-cliques straight from the index
    /tmp/nsi_base sample "$IDXF" $NQ >/tmp/e4q_${G}_${R}.txt 2>/dev/null
    NQGOT=$(wc -l </tmp/e4q_${G}_${R}.txt)
    # ours
    OURS=$(/tmp/nsi_query "$IDXF" bench /tmp/e4q_${G}_${R}.txt 2>/dev/null)
    OPT=$(echo "$OURS" | grep -oP 'point\(s=\d+\)=\K[0-9]+')
    OSP=$(echo "$OURS" | grep -oP 'spectrum\([0-9]+ cells\)=\K[0-9]+')
    # archive
    ARCH=$(/tmp/nsi_base bench "$IDXF" /tmp/e4q_${G}_${R}.txt $CAP 2>/dev/null)
    IMB=$(echo "$ARCH" | grep -oP '^index    : \K[0-9.]+(?= MB)' | head -1)
    IBR=$(echo "$ARCH" | grep -oP 'B / r-clique\)' -B0 >/dev/null; echo "$ARCH" | sed -n '2p' | grep -oP '\(\K[0-9.]+(?= B / r-clique)')
    AMB=$(echo "$ARCH" | grep -oP '^ARCHIVE  : \K[0-9.]+(?= MB)' | head -1)
    ABR=$(echo "$ARCH" | grep -oP 'ARCHIVE.*\(\K[0-9.]+(?= B / r-clique)' | head -1)
    SX=$(echo "$ARCH" | grep -oP 'INDEX IS \K[0-9.]+(?=x SMALLER)')
    APROBE=$(echo "$ARCH" | grep -oP 'ns/probe' -B0 >/dev/null; echo "$ARCH" | grep -oP '\K[0-9]+(?= ns/probe)')
    INFEAS=$(echo "$ARCH" | grep -c 'CANNOT BE MATERIALIZED')
    printf '%-20s %-4s %-10s %10s %12s %10s %12s %10s %10s %10s\n' \
      "$G" "$R" "$S0..$SMAX" "$BUILD" "${IMB:--}" "${IBR:--}" "${AMB:--}" "${ABR:--}" "${SX:--}" "${OPT:--}"
    {
      echo "### $G r=$R  (build ${BUILD}s, ${NQGOT} queries)"
      echo "$ARCH" | sed 's/^/    /'
      echo "    OURS  : point=${OPT:--} ns   spectrum=${OSP:--} ns"
      [ "$INFEAS" -gt 0 ] && echo "    VERDICT: archive infeasible -> index is the only sub-microsecond option"
      [ -n "${APROBE:-}" ] && echo "    QUERY RATIO: archive ${APROBE} ns / ours ${OPT} ns"
    } >>/data/wenqianz/e4_detail.out
  done
done
echo ""
echo "detail: /data/wenqianz/e4_detail.out"
echo "=== E4 INDEX DONE $(date) ==="
