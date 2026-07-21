#!/bin/bash
# §208 phase 2 on tods2:
#   PART A (G4) W = (Sum_P hostSz(P))/#r-cliques over every graph x cell we can reach.
#             Cheap: SCT_W_ONLY exits right after the maps phase, before supInit.
#   PART B    sparse footprints (SCT_SPARSE_FP) vs dense, TIME and MEMORY, 3 trials each.
#             3 trials because the local run showed sparse costing ~3% on ca-AstroPh and
#             saving ~2% on ca-HepPh: single runs cannot tell that from noise.
# Both parts matter jointly: the standing requirement is that speed AND memory both improve,
# so a config that wins one and loses the other is NOT acceptable and must be reported as such.
OUT=/data/wenqianz/sparse_w_tods2.out
exec >"$OUT" 2>&1
set -u
echo "=== SPARSE+W START $(date) ==="

REPO=/home/wenqianz/UNSW/pivoter
GD=/data/wenqianz
WANT_HEAD="${1:-}"

cd "$REPO" || { echo "FATAL: no repo"; exit 1; }
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then
  echo "FATAL: HEAD mismatch (want $WANT_HEAD)"; exit 1; fi

for sym in sparseFP leafFPSp spEqVec SCT_W_ONLY SCT_SPARSE_FP; do
  n=$(grep -c "$sym" region_native/region_native_sct_peel.cpp)
  echo "guard $sym = $n"
  [ "$n" -ge 1 ] || { echo "FATAL: $sym missing"; exit 1; }
done

cd "$REPO/region_native" || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn_sp region_native_sct_peel.cpp
BE=$?; echo "build_exit=$BE"; [ $BE -eq 0 ] || exit 1

export SCT_MAX_INC=200000000
export OMP_NUM_THREADS=1

echo ""
echo "=========================== PART A: W over all graphs ==========================="
printf '%-26s %-5s %12s %14s %12s %10s %8s\n' GRAPH CELL PATTERNS R_CLIQUES INCIDENCES COMPRESS W
echo "--------------------------------------------------------------------------------------------"
for G in ca-GrQc email-Eu-core ca-CondMat ca-HepPh ca-AstroPh com-amazon.ungraph com-dblp \
         web-Stanford web-Google web-it-2004 com-youtube wiki-Talk soc-pokec-relationships \
         com-dblp.p40 com-youtube.p40 web-Google.p40 web-it-2004.p40; do
  GR="$GD/$G.edges"; [ -f "$GR" ] || continue
  for CELL in "3 4" "3 5" "4 6" "5 7"; do
    set -- $CELL; R=$1; S=$2
    line=$(timeout 1200 env SCT_W_ONLY=1 /tmp/rn_sp "$GR" "$R" "$S" 2>&1 >/dev/null | grep '\[W ' )
    if [ -z "$line" ]; then printf '%-26s %-5s %s\n' "$G" "$R,$S" "(no result: timeout/guard/infeasible)"; continue; fi
    pa=$(echo "$line" | grep -oP 'patterns=\K[0-9]+')
    rc=$(echo "$line" | grep -oP 'r-cliques=\K[0-9]+')
    ic=$(echo "$line" | grep -oP 'incidences=\K[0-9]+')
    cp=$(echo "$line" | grep -oP 'compression=\K[0-9.]+')
    w=$(echo  "$line" | grep -oP 'W=\K[-0-9.]+')
    printf '%-26s %-5s %12s %14s %12s %10s %8s\n' "$G" "$R,$S" "$pa" "$rc" "$ic" "${cp}x" "$w"
  done
done

echo ""
echo "================ PART B: sparse vs dense, 3 trials, time AND memory ================"
printf '%-24s %-5s %-14s %10s %10s %8s %10s %10s %8s %s\n' \
  GRAPH CELL CONFIG T1 MEDIAN SPEEDUP RSS_MB MEM_SAVE BOTH CORE
echo "------------------------------------------------------------------------------------------------------------"
med3() { printf '%s\n' "$1" "$2" "$3" | sort -n | sed -n 2p; }

for spec in "ca-HepPh 3 5" "ca-AstroPh 4 6" "com-dblp 4 6" "web-Google 3 5" "com-youtube 3 5" "web-it-2004 3 4"; do
  set -- $spec; G=$1; R=$2; S=$3
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  declare -A TT; declare -A RR
  for cfg in dense sparse deconv both; do
    case $cfg in
      dense)  ENV="";;
      sparse) ENV="SCT_SPARSE_FP=1";;
      deconv) ENV="SCT_DECONV=1";;
      both)   ENV="SCT_SPARSE_FP=1 SCT_DECONV=1";;
    esac
    ts=(); ok=1
    for trial in 1 2 3; do
      timeout 2400 /usr/bin/time -v env $ENV /tmp/rn_sp "$GR" "$R" "$S" \
          >/tmp/w_${G}_${cfg}.out 2>/tmp/w_${G}_${cfg}.err || { ok=0; break; }
      t=$(grep -oP 'total=\K[0-9.]+' /tmp/w_${G}_${cfg}.out | tail -1)
      ts+=("$t")
    done
    [ $ok -eq 1 ] || { echo "$G ($R,$S) $cfg FAILED/TIMEOUT"; continue; }
    TT[$cfg]=$(med3 "${ts[0]}" "${ts[1]}" "${ts[2]}")
    RR[$cfg]=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/w_${G}_${cfg}.err)
    grep -E '^core=' /tmp/w_${G}_${cfg}.out >/tmp/w_${G}_${cfg}.core
    if [ "$cfg" = dense ]; then CORE="ref"; nb=$(wc -l </tmp/w_${G}_dense.core)
    elif [ "${nb:-0}" -eq 0 ]; then CORE="NO_REF"
    elif diff -q /tmp/w_${G}_dense.core /tmp/w_${G}_${cfg}.core >/dev/null; then CORE="BITEXACT"
    else CORE="**DIFFERS**"; fi
    sp=$(awk -v a="${TT[dense]:-0}" -v b="${TT[$cfg]}" 'BEGIN{if(b>0)printf "%.3f",a/b; else print "NA"}')
    ms=$(awk -v a="${RR[dense]:-0}" -v b="${RR[$cfg]}" 'BEGIN{if(b>0)printf "%.2f",a/b; else print "NA"}')
    both=$(awk -v s="$sp" -v m="$ms" 'BEGIN{ if(s>=1.0 && m>=1.0) print "BOTH-WIN"; else if(s<1.0 && m<1.0) print "BOTH-LOSE"; else print "SPLIT"; }')
    printf '%-24s %-5s %-14s %10s %10s %8s %10d %10s %8s %s\n' \
      "$G" "$R,$S" "$cfg" "${ts[0]}" "${TT[$cfg]}" "$sp" "$(( ${RR[$cfg]} / 1024 ))" "$ms" "$both" "$CORE"
  done
  unset TT RR nb
done

echo ""
echo "=== SPARSE+W DONE $(date) ==="
