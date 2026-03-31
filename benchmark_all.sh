#!/bin/bash
# =============================================================
# Full Benchmark: R=1..5, s from r+1 to max_clique
# 4 Ours + 3 REF algorithms, 4 parallel threads
# Timeout: 2 hours (7200 seconds)
# Output: benchmark_all_results.csv
# =============================================================

cd "$(dirname "$0")"

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
BIN="./build/bin/degeneracy_cliques"
TIMEOUT=7200
OUTCSV="benchmark_all_results.csv"
LOGDIR="bench_logs"
DATADIR="/data/wenqianz"
LOCKFILE="/tmp/bench_csv.lock"
NPROC=4

# Symlink graphs if not present
mkdir -p graphs
for f in com-dblp.edges com-youtube.edges web-Stanford.edges web-Google.edges soc-pokec-relationships.edges; do
  if [ ! -f "graphs/$f" ] && [ -f "$DATADIR/$f" ]; then
    ln -sf "$DATADIR/$f" "graphs/$f"
    echo "Linked graphs/$f -> $DATADIR/$f"
  fi
done

# Preserve existing results
if [ -f "$OUTCSV" ]; then
  cp "$OUTCSV" "${OUTCSV}.bak"
  echo "Backed up existing results to ${OUTCSV}.bak"
fi

# Build header + existing data
if [ ! -f "$OUTCSV" ] || ! head -1 "$OUTCSV" | grep -q 'graph,r,s'; then
  echo "graph,r,s,algorithm,time_ms,memory_kB,status" > "$OUTCSV"
fi
mkdir -p "$LOGDIR"

# run_one: runs a single (graph, r, s, algo) config
# Thread-safe: uses flock for CSV append, skips if result already exists
run_one() {
  local graph=$1 r=$2 s=$3 algo_label=$4 env_var=$5
  local logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_label}.log"

  # Skip if already has a result (OK or TIMEOUT)
  if grep -q "^${graph},${r},${s},${algo_label}," "$OUTCSV" 2>/dev/null; then
    return
  fi

  local result=""
  local exit_code=0
  result=$(env "${env_var}=1" timeout ${TIMEOUT}s $BIN "graphs/${graph}.edges" "$r" "$s" 2>&1) || exit_code=$?

  echo "$result" > "$logfile"

  local took=""
  local status="OK"
  if [ $exit_code -eq 124 ]; then
    took="TIMEOUT"
    status="TIMEOUT"
  else
    took=$(echo "$result" | grep -oP 'took: \K[0-9.]+(?= ms)' | tail -1) || true
    if [ -z "$took" ]; then
      took=$(echo "$result" | grep "^time:" | tail -1 | awk '{print $2}') || true
    fi
    if [ -z "$took" ]; then
      took="ERROR"
      status="ERROR(exit=$exit_code)"
    fi
  fi

  local mem=""
  mem=$(echo "$result" | grep -oP 'Final Memory: \s*\K[0-9]+') || true
  if [ -z "$mem" ]; then mem="N/A"; fi

  # Thread-safe CSV append
  (
    flock -x 200
    echo "$graph,$r,$s,$algo_label,$took,$mem,$status" >> "$OUTCSV"
  ) 200>"$LOCKFILE"

  printf "    %-12s r=%s s=%-3s %10sms (%s) mem=%s kB\n" "$algo_label" "$r" "$s" "$took" "$status" "$mem"
}
export -f run_one
export BIN TIMEOUT OUTCSV LOGDIR LOCKFILE

# ============ Step 2: Generate all jobs ============
echo "============================================================="
echo "  Full Benchmark (timeout=${TIMEOUT}s = 2h, ${NPROC} threads)"
echo "============================================================="

GRAPHS=(com-dblp com-youtube web-Stanford)
MAX_CLIQUE=(110 17 61)

JOBFILE=$(mktemp /tmp/bench_jobs.XXXXXX)

VALID_GRAPHS=()
VALID_MAX=()
for gi in "${!GRAPHS[@]}"; do
  if [ -f "graphs/${GRAPHS[$gi]}.edges" ]; then
    VALID_GRAPHS+=("${GRAPHS[$gi]}")
    VALID_MAX+=("${MAX_CLIQUE[$gi]}")
  else
    echo "SKIP: graphs/${GRAPHS[$gi]}.edges not found"
  fi
done

for gi in "${!VALID_GRAPHS[@]}"; do
  graph="${VALID_GRAPHS[$gi]}"
  omega="${VALID_MAX[$gi]}"

  for r in 1 2 3 4 5; do
    s_start=$((r + 1))
    if [ $s_start -gt $omega ]; then continue; fi

    for s in $(seq $s_start 1 $omega); do
      if [ $r -eq 1 ]; then
        echo "$graph $r $s Ours_ST PIVOTER_RUN_ST" >> "$JOBFILE"
        echo "$graph $r $s REF_R1 PIVOTER_NOOPT" >> "$JOBFILE"
      elif [ $r -eq 2 ]; then
        echo "$graph $r $s Ours_DCLP PIVOTER_RUN_R2_DCLP" >> "$JOBFILE"
        echo "$graph $r $s REF_R2 PIVOTER_NOOPT" >> "$JOBFILE"
      else
        echo "$graph $r $s Ours_V18 PIVOTER_RUN_ST_V18" >> "$JOBFILE"
        echo "$graph $r $s Ours_V19 PIVOTER_RUN_ST_V19" >> "$JOBFILE"
        echo "$graph $r $s Ours_V20 PIVOTER_RUN_ST_V20" >> "$JOBFILE"
        echo "$graph $r $s Ours_V11 PIVOTER_RUN_ST_V11" >> "$JOBFILE"
        echo "$graph $r $s REF_R3 PIVOTER_RUN_REF" >> "$JOBFILE"
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
cat "$JOBFILE" | xargs -P $NPROC -L 1 bash -c 'run_one $0 $1 $2 $3 $4'

rm -f "$JOBFILE" "$LOCKFILE"

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV ($(wc -l < "$OUTCSV") rows)"
echo "  Logs: $LOGDIR/"
echo "============================================================="
