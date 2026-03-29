#!/bin/bash
# =============================================================
# Full Benchmark: R=1,2,3,4 × multiple s × {Ours, Reference}
# Timeout: 20 minutes (1200 seconds)
# Output: bench_full_results.csv
# =============================================================

BIN="./build/bin/degeneracy_cliques"
TIMEOUT=1200
OUTCSV="bench_full_results.csv"
LOGDIR="bench_logs"
mkdir -p "$LOGDIR"

echo "graph,r,s,algorithm,time_ms,status" > "$OUTCSV"

run_one() {
  local graph=$1 r=$2 s=$3 algo_label=$4 env_var=$5
  local logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_label}.log"

  # Run with timeout
  local result
  result=$(eval "$env_var=1 perl -e 'alarm $TIMEOUT; exec @ARGV' $BIN graphs/${graph}.edges $r $s 2>&1")
  echo "$result" > "$logfile"

  # Extract time
  local took=$(echo "$result" | grep -E "took:.*ms" | tail -1 | sed 's/.*took: //' | sed 's/ ms//')
  local status="OK"

  if [ -z "$took" ]; then
    # Try "time:" fallback
    took=$(echo "$result" | grep "^time:" | tail -1 | awk '{print $2}')
  fi

  if [ -z "$took" ]; then
    took="TIMEOUT"
    status="TIMEOUT"
  fi

  echo "$graph,$r,$s,$algo_label,$took,$status" >> "$OUTCSV"
  echo "  $algo_label: ${took}ms ($status)"
}

# =============================================================
# Graph configs: graph_name, max_clique_size
# =============================================================
declare -a GRAPHS=("email-Eu-core" "com-dblp" "web-Stanford" "com-youtube")
declare -a MAX_CLIQUE=(18 110 60 17)

echo "============================================================="
echo "  FULL BENCHMARK (timeout=${TIMEOUT}s)"
echo "============================================================="
echo ""

for gi in "${!GRAPHS[@]}"; do
  graph="${GRAPHS[$gi]}"
  omega="${MAX_CLIQUE[$gi]}"
  echo "=== $graph (ω=$omega) ==="

  # ----------------------------------------------------------
  # R=1: Ours (ST_V2) vs Reference
  # s from 3 to min(omega, 30)
  # ----------------------------------------------------------
  echo "--- R=1 ---"
  max_s=$((omega < 30 ? omega : 30))
  for s in $(seq 3 2 $max_s); do
    echo -n "  r=1 s=$s: "
    run_one "$graph" 1 "$s" "R1_ST" "PIVOTER_RUN_ST"
    run_one "$graph" 1 "$s" "R1_REF" "PIVOTER_RUN_REF"
  done

  # ----------------------------------------------------------
  # R=2: Ours (DCLP) vs Reference
  # s from 3 to min(omega, 30)
  # ----------------------------------------------------------
  echo "--- R=2 ---"
  for s in $(seq 3 2 $max_s); do
    echo -n "  r=2 s=$s: "
    run_one "$graph" 2 "$s" "R2_DCLP" "PIVOTER_RUN_R2_DCLP"
    run_one "$graph" 2 "$s" "R2_REF" "PIVOTER_RUN_REF"
  done

  # ----------------------------------------------------------
  # R=3: Ours (V18) vs Reference
  # s from 5 to min(omega, 30)
  # ----------------------------------------------------------
  echo "--- R=3 ---"
  max_s3=$((omega < 30 ? omega : 30))
  for s in $(seq 5 2 $max_s3); do
    echo -n "  r=3 s=$s: "
    run_one "$graph" 3 "$s" "R3_V18" "PIVOTER_RUN_ST_V18"
    run_one "$graph" 3 "$s" "R3_V11" "PIVOTER_RUN_ST_V11"
    run_one "$graph" 3 "$s" "R3_REF" "PIVOTER_RUN_REF"
  done

  # ----------------------------------------------------------
  # R=4: Ours (V18) vs Reference
  # s from 6 to min(omega, 20)
  # ----------------------------------------------------------
  echo "--- R=4 ---"
  max_s4=$((omega < 20 ? omega : 20))
  for s in $(seq 6 2 $max_s4); do
    echo -n "  r=4 s=$s: "
    run_one "$graph" 4 "$s" "R4_V18" "PIVOTER_RUN_ST_V18"
    run_one "$graph" 4 "$s" "R4_REF" "PIVOTER_RUN_REF"
  done

  echo ""
done

echo "============================================================="
echo "  DONE. Results in $OUTCSV"
echo "============================================================="
