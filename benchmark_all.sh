#!/bin/bash
# =============================================================
# Full Benchmark: R=1..5, s from r+1 to max_clique
# 4 Ours + 3 REF algorithms, 4 parallel threads
# Timeout: 1 hour (3600 seconds)
# Output: benchmark_all_results.csv
#
# Usage:
#   bash benchmark_all.sh          # normal: generate jobs & run parallel
#   bash benchmark_all.sh --run G R S LABEL ENV   # internal: run single job
# =============================================================

cd "$(dirname "$0")"

BIN="./build/bin/degeneracy_cliques"
TIMEOUT=3600
OUTCSV="benchmark_all_results.csv"
LOGDIR="bench_logs"
DATADIR="/data/wenqianz"
LOCKFILE="/tmp/bench_csv.lock"
NPROC=8

# ============ Internal: run a single job ============
if [ "$1" = "--run" ]; then
  graph=$2; r=$3; s=$4; algo_label=$5; env_var=$6
  logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_label}.log"

  # Skip if already has a result
  if grep -q "^${graph},${r},${s},${algo_label}," "$OUTCSV" 2>/dev/null; then
    exit 0
  fi

  result=""
  exit_code=0
  result=$(env "${env_var}=1" timeout ${TIMEOUT}s $BIN "graphs/${graph}.edges" "$r" "$s" 2>&1) || exit_code=$?
  echo "$result" > "$logfile"

  took=""
  status="OK"
  if [ $exit_code -eq 124 ]; then
    took="TIMEOUT"; status="TIMEOUT"
  else
    # Priority: NucleusCoreDecomposition took > ST_V* took > time: > any took:
    took=$(echo "$result" | grep -oP 'NucleusCoreDecomposition took: \K[0-9.]+(?= ms)') || true
    if [ -z "$took" ]; then
      took=$(echo "$result" | grep -oP 'ST_V\S+ took: \K[0-9.]+(?= ms)') || true
    fi
    if [ -z "$took" ]; then
      took=$(echo "$result" | grep "^time:" | tail -1 | awk '{print $2}') || true
    fi
    if [ -z "$took" ]; then
      took="ERROR"; status="ERROR(exit=$exit_code)"
    fi
  fi

  mem=$(echo "$result" | grep -oP 'Final Memory: \s*\K[0-9]+') || true
  if [ -z "$mem" ]; then mem="N/A"; fi

  (
    flock -x 200
    echo "$graph,$r,$s,$algo_label,$took,$mem,$status" >> "$OUTCSV"
  ) 200>"$LOCKFILE"

  printf "    %-12s r=%s s=%-3s %10sms (%s) mem=%s kB\n" "$algo_label" "$r" "$s" "$took" "$status" "$mem"
  exit 0
fi

# ============ Step 0: Build ============
echo "============================================================="
echo "  Step 0: Clean build"
echo "============================================================="
rm -rf build
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build -j 12 --target degeneracy_cliques 2>&1 | tail -5
echo "Build done."
echo ""

# ============ Step 1: Setup ============
mkdir -p graphs
for f in com-dblp.edges com-youtube.edges web-Stanford.edges web-Google.edges soc-pokec-relationships.edges; do
  if [ ! -f "graphs/$f" ] && [ -f "$DATADIR/$f" ]; then
    ln -sf "$DATADIR/$f" "graphs/$f"
    echo "Linked graphs/$f -> $DATADIR/$f"
  fi
done

if [ -f "$OUTCSV" ]; then
  cp "$OUTCSV" "${OUTCSV}.bak"
  echo "Backed up existing results to ${OUTCSV}.bak"
fi
if [ ! -f "$OUTCSV" ] || ! head -1 "$OUTCSV" | grep -q 'graph,r,s'; then
  echo "graph,r,s,algorithm,time_ms,memory_kB,status" > "$OUTCSV"
fi
mkdir -p "$LOGDIR"

# ============ Step 2: Generate all jobs ============
echo "============================================================="
echo "  Full Benchmark (timeout=${TIMEOUT}s = 1h, ${NPROC} threads)"
echo "============================================================="

GRAPHS=(com-dblp com-youtube web-Stanford)
MAX_CLIQUE=(110 17 61)

JOBFILE=$(mktemp /tmp/bench_jobs.XXXXXX)

for gi in "${!GRAPHS[@]}"; do
  graph="${GRAPHS[$gi]}"
  omega="${MAX_CLIQUE[$gi]}"
  if [ ! -f "graphs/${graph}.edges" ]; then
    echo "SKIP: graphs/${graph}.edges not found"
    continue
  fi

  for r in 1 2 3 4 5; do
    s_start=$((r + 1))
    if [ $s_start -gt $omega ]; then continue; fi

    for s in $(seq $s_start 1 $omega); do
      if [ $r -eq 1 ]; then
        echo "--run $graph $r $s Ours_ST PIVOTER_RUN_ST" >> "$JOBFILE"
        echo "--run $graph $r $s REF_R1 PIVOTER_NOOPT" >> "$JOBFILE"
      elif [ $r -eq 2 ]; then
        echo "--run $graph $r $s Ours_DCLP PIVOTER_RUN_R2_DCLP" >> "$JOBFILE"
        echo "--run $graph $r $s REF_R2 PIVOTER_NOOPT" >> "$JOBFILE"
      else
        echo "--run $graph $r $s Ours_V18 PIVOTER_RUN_ST_V18" >> "$JOBFILE"
        echo "--run $graph $r $s Ours_V20 PIVOTER_RUN_ST_V20" >> "$JOBFILE"
        echo "--run $graph $r $s REF_R3 PIVOTER_RUN_REF" >> "$JOBFILE"
      fi
    done
  done
done

TOTAL=$(wc -l < "$JOBFILE")
EXISTING=$(tail -n +2 "$OUTCSV" 2>/dev/null | wc -l | tr -d ' ')
echo "Total jobs: $TOTAL, existing results: $EXISTING"
echo "Running with $NPROC parallel threads..."
echo ""

# ============ Step 3: Run in parallel ============
cat "$JOBFILE" | xargs -P $NPROC -L 1 bash "$(realpath "$0")"

rm -f "$JOBFILE" "$LOCKFILE"

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV ($(wc -l < "$OUTCSV") rows)"
echo "  Logs: $LOGDIR/"
echo "============================================================="
