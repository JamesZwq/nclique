#!/bin/bash
# =============================================================
# V3 Phase Benchmark: ST vs V2 vs V3 with phase breakdown
# Runs on tods2 server. 8 parallel threads. 10 min timeout.
#
# Usage:
#   bash bench_v3_phases.sh              # generate jobs + run
#   bash bench_v3_phases.sh --run ...    # internal: single job
# =============================================================

cd "$(dirname "$0")"

BIN="./build/bin/degeneracy_cliques"
TIMEOUT=600
OUTCSV="bench_v3_phases_results.csv"
LOGDIR="bench_v3_phases_logs"
DATADIR="/data/wenqianz"
LOCKFILE="/tmp/bench_v3_phases.lock"
NPROC=8

# ============ Internal: run a single job ============
if [ "$1" = "--run" ]; then
  graph=$2; r=$3; s=$4; algo=$5; env_var=$6
  logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo}.log"

  # Skip if done
  if grep -q "^${graph},${r},${s},${algo}," "$OUTCSV" 2>/dev/null; then
    exit 0
  fi

  result=""
  exit_code=0
  result=$(env "${env_var}=1" timeout ${TIMEOUT}s $BIN "graphs/${graph}.edges" "$r" "$s" 2>&1) || exit_code=$?
  echo "$result" > "$logfile"

  status="OK"
  if [ $exit_code -eq 124 ]; then
    status="TIMEOUT"
  elif [ $exit_code -ne 0 ]; then
    status="ERROR($exit_code)"
  fi

  # Extract timing
  total=$(echo "$result" | grep -oP 'NucleusCoreDecomposition took: \K[0-9.]+(?= ms)') || true

  # Extract memory
  graph_mem=$(echo "$result" | grep "Graph Memory" | grep -oP '[0-9]+' | tail -1) || true
  index_mem=$(echo "$result" | grep "Other Index Memory\|Other index" | grep -oP '[0-9]+' | tail -1) || true
  final_mem=$(echo "$result" | grep "Final Memory" | grep -oP '[0-9]+' | tail -1) || true

  # V3-specific phases
  sdct_time=$(echo "$result" | grep -oP 'SDCT_MaxClique took: \K[0-9.]+(?= ms)') || true
  maxcliq_time=$(echo "$result" | grep -oP 'MaxCliqEnum \(V3\) took: \K[0-9.]+(?= ms)') || true
  merge_time=$(echo "$result" | grep -oP 'r-Mergeable classification: \K[0-9.]+(?= ms)') || true
  cpi_time=$(echo "$result" | grep "CPI counting time:" | sed 's/.*: //' | sed 's/ ms//' | tr -d ' ') || true
  pathinfo_time=$(echo "$result" | grep "PathInfo build time:" | sed 's/.*: //' | sed 's/ ms//' | tr -d ' ') || true
  peel_time=$(echo "$result" | grep "Peeling time:" | sed 's/.*: //' | sed 's/ ms//' | tr -d ' ') || true
  v3_total=$(echo "$result" | grep "Total time:" | sed 's/.*: //' | sed 's/ ms//' | tr -d ' ') || true

  # V3-specific stats
  fully_merge=$(echo "$result" | grep "Fully mergeable" | grep -oP '\d+(?= \()') || true
  remaining=$(echo "$result" | grep "Remaining regions:" | sed 's/.*: //' | tr -d ' ') || true
  classes=$(echo "$result" | grep "Overlap classes:" | sed 's/.*: //' | tr -d ' ') || true
  tuples=$(echo "$result" | grep "r-tuples:" | sed 's/.*: //' | tr -d ' ') || true
  recursive=$(echo "$result" | grep "Total recursive calls:" | sed 's/.*: //' | tr -d ' ') || true
  max_core=$(echo "$result" | grep "Max core:" | sed 's/.*: //' | head -1 | tr -d ' ') || true

  # V2-specific
  v2_enum=$(echo "$result" | grep "Enumeration time:" | sed 's/.*: //' | sed 's/ ms//' | tr -d ' ') || true
  v2_stuples=$(echo "$result" | grep "s-tuples:" | sed 's/.*: //' | head -1 | tr -d ' ') || true
  v2_peel=$(echo "$result" | grep "Peeling time:" | sed 's/.*: //' | sed 's/ ms//' | tr -d ' ') || true

  # ST-specific
  st_fused=$(echo "$result" | grep -oP 'fused build\+counting.*took: \K[0-9.]+(?= ms)') || true
  st_bk=$(echo "$result" | grep "BK:" | sed 's/.*BK: *//' | tr -d ' ') || true

  (
    flock -x 200
    echo "${graph},${r},${s},${algo},${status},${total:-NA},${graph_mem:-NA},${index_mem:-NA},${final_mem:-NA},${sdct_time:-NA},${maxcliq_time:-NA},${merge_time:-NA},${cpi_time:-NA},${pathinfo_time:-NA},${peel_time:-NA},${v3_total:-NA},${fully_merge:-NA},${remaining:-NA},${classes:-NA},${tuples:-NA},${recursive:-NA},${max_core:-NA},${v2_enum:-NA},${v2_stuples:-NA},${v2_peel:-NA},${st_fused:-NA},${st_bk:-NA}" >> "$OUTCSV"
  ) 200>"$LOCKFILE"

  printf "  %-12s %-18s r=%s s=%-2s total=%sms peel=%sms status=%s\n" \
    "$algo" "$graph" "$r" "$s" "${total:-NA}" "${peel_time:-NA}" "$status"
  exit 0
fi

# ============ Step 0: Build ============
echo "============================================================="
echo "  V3 Phase Benchmark"
echo "============================================================="
echo ""

# Link data from server data dir
mkdir -p graphs
for f in com-dblp.edges com-youtube.edges web-Stanford.edges web-Google.edges \
         soc-pokec-relationships.edges web-it-2004.edges \
         dblp-core30.edges email-Eu-core.edges; do
  if [ ! -f "graphs/$f" ] && [ -f "$DATADIR/$f" ]; then
    ln -sf "$DATADIR/$f" "graphs/$f"
    echo "Linked $f"
  fi
done

# Build
echo "Building..."
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -1
cmake --build build -j 12 --target degeneracy_cliques 2>&1 | tail -1
echo "Build done."
echo ""

# Setup output
if [ ! -f "$OUTCSV" ] || ! head -1 "$OUTCSV" | grep -q 'graph,r,s'; then
  echo "graph,r,s,algo,status,total_ms,graph_mem_kB,index_mem_kB,final_mem_kB,sdct_ms,maxcliq_ms,merge_ms,cpi_ms,pathinfo_ms,peel_ms,v3_total_ms,fully_merge,remaining,classes,tuples,recursive_calls,max_core,v2_enum_ms,v2_stuples,v2_peel_ms,st_fused_ms,st_bk_ms" > "$OUTCSV"
fi
mkdir -p "$LOGDIR"

# ============ Generate jobs ============
GRAPHS=(com-dblp web-Stanford web-it-2004 dblp-core30 email-Eu-core)
ALGOS="ST V2 V3"

JOBFILE=$(mktemp /tmp/bench_v3_jobs.XXXXXX)

for graph in "${GRAPHS[@]}"; do
  if [ ! -f "graphs/${graph}.edges" ]; then
    echo "SKIP: graphs/${graph}.edges not found"
    continue
  fi

  for rs in "3 4" "3 5" "3 6" "3 8" "3 10" "4 5" "4 6" "5 6"; do
    r=${rs% *}; s=${rs#* }

    for algo in $ALGOS; do
      case $algo in
        ST) env_var="PIVOTER_RUN_ST" ;;
        V2) env_var="PIVOTER_RUN_REGION_V2" ;;
        V3) env_var="PIVOTER_RUN_REGION_V3" ;;
      esac
      echo "--run $graph $r $s $algo $env_var" >> "$JOBFILE"
    done
  done
done

TOTAL=$(wc -l < "$JOBFILE")
EXISTING=$(tail -n +2 "$OUTCSV" 2>/dev/null | wc -l | tr -d ' ')
echo "Total jobs: $TOTAL, existing: $EXISTING, parallel: $NPROC"
echo "Timeout per job: ${TIMEOUT}s"
echo "Output: $OUTCSV"
echo ""

# ============ Run ============
cat "$JOBFILE" | xargs -P $NPROC -L 1 bash "$(realpath "$0")"

rm -f "$JOBFILE" "$LOCKFILE"

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV ($(wc -l < "$OUTCSV") rows)"
echo "  Logs: $LOGDIR/"
echo "============================================================="

# ============ Summary table ============
echo ""
echo "=== Summary ==="
echo ""
printf "%-18s %3s %3s | %8s %8s %8s | %8s %8s %8s %8s\n" \
  "Graph" "r" "s" "ST" "V2" "V3" "V3:CPI" "V3:peel" "V3:merge" "V3:tuples"
echo "------------------------------------------------------------------------------------"

tail -n +2 "$OUTCSV" | sort -t',' -k1,1 -k2,2n -k3,3n | awk -F',' '
{
  key=$1","$2","$3
  if ($4=="ST") st[key]=$6
  if ($4=="V2") v2[key]=$6
  if ($4=="V3") {
    v3[key]=$6; v3cpi[key]=$13; v3peel[key]=$15; v3merge[key]=$17; v3tuples[key]=$20
  }
  seen[key]=1
}
END {
  for (key in seen) {
    split(key, a, ",")
    printf "%-18s %3s %3s | %8s %8s %8s | %8s %8s %8s %8s\n",
      a[1], a[2], a[3],
      (st[key]?st[key]:"NA"), (v2[key]?v2[key]:"NA"), (v3[key]?v3[key]:"NA"),
      (v3cpi[key]?v3cpi[key]:"NA"), (v3peel[key]?v3peel[key]:"NA"),
      (v3merge[key]?v3merge[key]:"NA"), (v3tuples[key]?v3tuples[key]:"NA")
  }
}' | sort
