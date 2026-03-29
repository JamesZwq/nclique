#!/bin/bash
# =============================================================
# Full Benchmark: R=1,2,3,4 × multiple s × {Ours, Reference}
# Run on server: ssh tods2, data at /data/wenqianz/
# Timeout: 20 minutes (1200 seconds)
# Output: bench_full_results.csv
# =============================================================

set -e

# ============ Step 0: Build ============
echo "============================================================="
echo "  Step 0: Clean build"
echo "============================================================="
cd "$(dirname "$0")"
rm -rf build
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build -j 12 --target degeneracy_cliques 2>&1 | tail -5
echo "Build done."
echo ""

# ============ Step 1: Setup ============
BIN="./build/bin/degeneracy_cliques"
TIMEOUT=1200
OUTCSV="bench_full_results.csv"
LOGDIR="bench_logs"
DATADIR="/data/wenqianz"

# Symlink graphs if not present
mkdir -p graphs
for f in com-dblp.edges com-youtube.edges web-Stanford.edges web-Google.edges soc-pokec-relationships.edges; do
  if [ ! -f "graphs/$f" ] && [ -f "$DATADIR/$f" ]; then
    ln -sf "$DATADIR/$f" "graphs/$f"
    echo "Linked graphs/$f -> $DATADIR/$f"
  fi
done

echo "graph,r,s,algorithm,time_ms,status" > "$OUTCSV"
mkdir -p "$LOGDIR"

run_one() {
  local graph=$1 r=$2 s=$3 algo_label=$4 env_var=$5
  local logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_label}.log"
  local result
  result=$(eval "$env_var=1 perl -e 'alarm $TIMEOUT; exec @ARGV' $BIN graphs/${graph}.edges $r $s 2>&1")
  echo "$result" > "$logfile"

  local took=$(echo "$result" | grep -E "took:.*ms" | tail -1 | sed 's/.*took: //' | sed 's/ ms//')
  local status="OK"
  if [ -z "$took" ]; then
    took=$(echo "$result" | grep "^time:" | tail -1 | awk '{print $2}')
  fi
  if [ -z "$took" ]; then
    took="TIMEOUT"
    status="TIMEOUT"
  fi

  echo "$graph,$r,$s,$algo_label,$took,$status" >> "$OUTCSV"
  printf "    %-10s %s (%s)\n" "$algo_label" "${took}ms" "$status"
}

# ============ Step 2: Run experiments ============
echo "============================================================="
echo "  Full Benchmark (timeout=${TIMEOUT}s = 20min)"
echo "============================================================="

# Graph configs: name, max_clique
GRAPHS=(com-dblp com-youtube web-Stanford web-Google)
MAX_CLIQUE=(110 17 60 40)

# Check which graphs exist
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
  echo ""
  echo "=== $graph (ω=$omega) ==="

  # R=1: s = 3,5,7,...,min(omega,30)
  echo "  --- R=1 ---"
  max_s=$((omega < 30 ? omega : 30))
  for s in $(seq 3 2 $max_s); do
    echo "  r=1 s=$s:"
    run_one "$graph" 1 "$s" "R1_Ours" "PIVOTER_RUN_ST"
    run_one "$graph" 1 "$s" "R1_REF" "PIVOTER_RUN_REF"
  done

  # R=2: s = 3,5,7,...,min(omega,30)
  echo "  --- R=2 ---"
  for s in $(seq 3 2 $max_s); do
    echo "  r=2 s=$s:"
    run_one "$graph" 2 "$s" "R2_DCLP" "PIVOTER_RUN_R2_DCLP"
    run_one "$graph" 2 "$s" "R2_REF" "PIVOTER_RUN_REF"
  done

  # R=3: s = 5,7,9,...,min(omega,30)
  echo "  --- R=3 ---"
  max_s3=$((omega < 30 ? omega : 30))
  for s in $(seq 5 2 $max_s3); do
    echo "  r=3 s=$s:"
    run_one "$graph" 3 "$s" "R3_V18" "PIVOTER_RUN_ST_V18"
    run_one "$graph" 3 "$s" "R3_V11" "PIVOTER_RUN_ST_V11"
    run_one "$graph" 3 "$s" "R3_REF" "PIVOTER_RUN_REF"
  done

  # R=4: s = 6,8,10,...,min(omega,20)
  echo "  --- R=4 ---"
  max_s4=$((omega < 20 ? omega : 20))
  for s in $(seq 6 2 $max_s4); do
    echo "  r=4 s=$s:"
    run_one "$graph" 4 "$s" "R4_V18" "PIVOTER_RUN_ST_V18"
    run_one "$graph" 4 "$s" "R4_REF" "PIVOTER_RUN_REF"
  done
done

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV"
echo "  Logs: $LOGDIR/"
echo "============================================================="
