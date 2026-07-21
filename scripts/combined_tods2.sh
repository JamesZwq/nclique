#!/bin/bash
# §210 phase 3 on tods2: the CLAMP PRE-FILTER across graphs, and the full combination.
#
# Two rules this script enforces, both learned the hard way tonight (§210):
#  1. ATTRIBUTE THE METRIC TO THE PHASE YOU CHANGED. All three optimisations live inside the peel.
#     Comparing totals on MCE-dominated graphs manufactured both a false win (com-youtube 1.027x)
#     and a false loss (web-it-2004 0.943x) from the same MCE noise. PEEL is the primary column;
#     TOTAL is reported beside it, and the untouched phases are printed so they can be checked equal.
#  2. THREE TRIALS, MEDIAN. Single runs cannot separate the 3-8% effects that decide these calls.
#
# Standing requirement: speed AND memory must BOTH improve. Rows are labelled accordingly.
OUT=/data/wenqianz/combined_tods2.out
exec >"$OUT" 2>&1
set -u
echo "=== COMBINED START $(date) ==="

REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz; WANT_HEAD="${1:-}"
cd "$REPO" || exit 1
git fetch -q origin && git reset -q --hard origin/main
HEAD_NOW=$(git log --oneline -1); echo "HEAD: $HEAD_NOW"
if [ -n "$WANT_HEAD" ] && ! echo "$HEAD_NOW" | grep -q "$WANT_HEAD"; then
  echo "FATAL: HEAD mismatch (want $WANT_HEAD)"; exit 1; fi
for sym in clampPF patLive SCT_CLAMP_PF SCT_SPARSE_FP SCT_DECONV; do
  n=$(grep -c "$sym" region_native/region_native_sct_peel.cpp); echo "guard $sym = $n"
  [ "$n" -ge 1 ] || { echo "FATAL: $sym missing"; exit 1; }; done
cd region_native || exit 1
g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn_cmb region_native_sct_peel.cpp
BE=$?; echo "build_exit=$BE"; [ $BE -eq 0 ] || exit 1

export SCT_MAX_INC=200000000 OMP_NUM_THREADS=1
med3(){ printf '%s\n' "$1" "$2" "$3" | sort -n | sed -n 2p; }

printf '%-22s %-5s %-16s %9s %9s %8s %9s %8s %9s %-9s %s\n' \
  GRAPH CELL CONFIG PEEL_MED TOT_MED PEEL_SPD RSS_MB MEM_SAVE MCE VERDICT CORE
echo "--------------------------------------------------------------------------------------------------------------------"

for spec in "ca-GrQc 3 5" "ca-CondMat 3 5" "email-Eu-core 3 5" "ca-HepPh 3 5" "ca-AstroPh 4 6" \
            "com-dblp 4 6" "com-amazon.ungraph 3 5" "web-Google 3 5" "com-youtube 3 5" "web-it-2004 3 4"; do
  set -- $spec; G=$1; R=$2; S=$3
  GR="$GD/$G.edges"; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  refpeel=""; refrss=""
  for cfg in dense pf sparse+deconv all; do
    case $cfg in
      dense)          ENV="";;
      pf)             ENV="SCT_CLAMP_PF=1";;
      sparse+deconv)  ENV="SCT_SPARSE_FP=1 SCT_DECONV=1";;
      all)            ENV="SCT_CLAMP_PF=1 SCT_SPARSE_FP=1 SCT_DECONV=1";;
    esac
    p1=""; p2=""; p3=""; t1=""; t2=""; t3=""; fail=0
    for i in 1 2 3; do
      timeout 3000 /usr/bin/time -v env $ENV /tmp/rn_cmb "$GR" "$R" "$S" \
        >/tmp/c_${G}_${cfg}.out 2>/tmp/c_${G}_${cfg}.err || { fail=1; break; }
      eval "p$i=\$(grep -oP 'peel=\K[0-9.]+'  /tmp/c_${G}_${cfg}.out | tail -1)"
      eval "t$i=\$(grep -oP 'total=\K[0-9.]+' /tmp/c_${G}_${cfg}.out | tail -1)"
    done
    [ $fail -eq 0 ] || { printf '%-22s %-5s %-16s %s\n' "$G" "$R,$S" "$cfg" "FAILED/TIMEOUT"; continue; }
    pm=$(med3 "$p1" "$p2" "$p3"); tm=$(med3 "$t1" "$t2" "$t3")
    rss=$(grep -oP 'Maximum resident set size \(kbytes\): \K[0-9]+' /tmp/c_${G}_${cfg}.err)
    mce=$(grep -oP 'MCE=\K[0-9.]+' /tmp/c_${G}_${cfg}.out | tail -1)
    grep -E '^core=' /tmp/c_${G}_${cfg}.out >/tmp/c_${G}_${cfg}.core
    if [ "$cfg" = dense ]; then refpeel=$pm; refrss=$rss; CORE=ref
    else
      nb=$(wc -l </tmp/c_${G}_dense.core)
      if [ "$nb" -eq 0 ]; then CORE=NO_REF
      elif diff -q /tmp/c_${G}_dense.core /tmp/c_${G}_${cfg}.core >/dev/null; then CORE=BITEXACT
      else CORE='**DIFFERS**'; fi
    fi
    sp=$(awk -v a="$refpeel" -v b="$pm"  'BEGIN{if(b>0)printf "%.3f",a/b; else print "-"}')
    ms=$(awk -v a="$refrss"  -v b="$rss" 'BEGIN{if(b>0)printf "%.2f",a/b; else print "-"}')
    vd=$(awk -v s="$sp" -v m="$ms" 'BEGIN{ if(s=="-"||m=="-") print "-"; else if(s>=0.99&&m>=0.99) print "BOTH-WIN"; else if(s<0.99&&m<0.99) print "BOTH-LOSE"; else print "SPLIT"; }')
    [ "$cfg" = dense ] && vd=ref
    printf '%-22s %-5s %-16s %9s %9s %8s %9d %8s %9s %-9s %s\n' \
      "$G" "$R,$S" "$cfg" "$pm" "$tm" "$sp" "$((rss/1024))" "$ms" "$mce" "$vd" "$CORE"
  done
  echo ""
done
echo "=== COMBINED DONE $(date) ==="
