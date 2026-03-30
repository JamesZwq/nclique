#!/bin/bash
# =============================================================
# Full Benchmark: R=1..5, s from r+1 to max_clique
# 4 Ours + 3 REF algorithms
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

# Symlink graphs if not present
mkdir -p graphs
for f in com-dblp.edges com-youtube.edges web-Stanford.edges web-Google.edges soc-pokec-relationships.edges; do
  if [ ! -f "graphs/$f" ] && [ -f "$DATADIR/$f" ]; then
    ln -sf "$DATADIR/$f" "graphs/$f"
    echo "Linked graphs/$f -> $DATADIR/$f"
  fi
done

echo "graph,r,s,algorithm,time_ms,memory_kB,status" > "$OUTCSV"
mkdir -p "$LOGDIR"

run_one() {
  local graph=$1 r=$2 s=$3 algo_label=$4 env_var=$5
  local logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_label}.log"
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

  echo "$graph,$r,$s,$algo_label,$took,$mem,$status" >> "$OUTCSV"
  printf "    %-12s %10sms (%s) mem=%s kB\n" "$algo_label" "$took" "$status" "$mem"
}

# ============ Step 2: Run experiments ============
echo "============================================================="
echo "  Full Benchmark (timeout=${TIMEOUT}s = 2h)"
echo "============================================================="

# Graph configs: name, max_clique (omega)
GRAPHS=(com-dblp com-youtube web-Stanford)
MAX_CLIQUE=(110 17 61)

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

  for r in 1 2 3 4 5; do
    # s ranges from r+1 to min(omega, cap)
    # Cap s to avoid extremely long runs
    local_cap=$omega
    s_start=$((r + 1))
    if [ $s_start -gt $local_cap ]; then continue; fi

    echo "  --- R=$r ---"
    for s in $(seq $s_start 1 $local_cap); do
      echo "  r=$r s=$s:"

      if [ $r -eq 1 ]; then
        run_one "$graph" "$r" "$s" "Ours_ST" "PIVOTER_RUN_ST"
        run_one "$graph" "$r" "$s" "REF_R1" "PIVOTER_NOOPT"
      elif [ $r -eq 2 ]; then
        run_one "$graph" "$r" "$s" "Ours_DCLP" "PIVOTER_RUN_R2_DCLP"
        run_one "$graph" "$r" "$s" "REF_R2" "PIVOTER_NOOPT"
      else
        # R>=3: run V18, V11, and REF
        run_one "$graph" "$r" "$s" "Ours_V18" "PIVOTER_RUN_ST_V18"
        run_one "$graph" "$r" "$s" "Ours_V11" "PIVOTER_RUN_ST_V11"
        run_one "$graph" "$r" "$s" "REF_R3" "PIVOTER_RUN_REF"
      fi
    done
  done
done

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV"
echo "  Logs: $LOGDIR/"
echo "============================================================="
