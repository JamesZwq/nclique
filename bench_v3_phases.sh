#!/bin/bash
# =============================================================
# V3 Phase Benchmark: ST vs V2 vs V3 with phase breakdown
# Runs on tods2 server. 8 parallel threads. 10 min timeout.
#
# Usage:
#   nohup bash bench_v3_phases.sh > bench_v3_phases.log 2>&1 &
#   bash bench_v3_phases.sh --run ...    # internal: single job
# =============================================================

cd "$(dirname "$0")"

BIN="./build/bin/degeneracy_cliques"
TIMEOUT=600
OUTCSV="bench_v3_phases_results.csv"
LOGDIR="bench_v3_phases_logs"
DATADIR="/data/wenqianz"
LOCKFILE="/tmp/bench_v3_phases.lock"
MAX_NPROC=32                 # max parallel jobs
MEM_LIMIT_GB=300             # don't launch new jobs if total usage > this
MEM_KILL_GB=450              # kill newest job if usage exceeds this
MEM_CHECK_INTERVAL=5         # seconds between checks when waiting
LAUNCH_SETTLE=3              # brief pause after launch to let process start

# ============ Internal: run a single job ============
if [ "$1" = "--run" ]; then
  graph=$2; r=$3; s=$4; algo=$5; env_var=$6
  logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo}.log"

  # Skip if done
  if grep -q "^${graph},${r},${s},${algo}," "$OUTCSV" 2>/dev/null; then
    exit 0
  fi

  # Trap signals: kill child process on TERM/INT (so parent kill propagates)
  CHILD_PID=""
  trap 'if [ -n "$CHILD_PID" ]; then kill "$CHILD_PID" 2>/dev/null; wait "$CHILD_PID" 2>/dev/null; fi; exit 143' TERM INT

  result=""
  exit_code=0
  extra_env=""
  if [ "$algo" = "V3_NP" ]; then
    extra_env="PIVOTER_V3_NO_PRIVATE=1"
  fi
  env "${env_var}=1" $extra_env timeout ${TIMEOUT}s $BIN "graphs/${graph}.edges" "$r" "$s" > "$logfile" 2>&1 &
  CHILD_PID=$!
  wait "$CHILD_PID"
  exit_code=$?
  CHILD_PID=""
  result=$(cat "$logfile")

  # Determine status
  status="OK"
  if [ $exit_code -eq 124 ]; then
    status="TIMEOUT"
  elif [ $exit_code -ne 0 ]; then
    status="ERROR($exit_code)"
  fi

  # Quick extract for console output only (Python does the real parsing)
  total=$(grep 'NucleusCoreDecomposition took:' "$logfile" 2>/dev/null | awk '{print $3}') || true
  peel=$(grep 'Peeling time:' "$logfile" 2>/dev/null | tail -1 | sed 's/.*Peeling time:[[:space:]]*//' | sed 's/[[:space:]]*ms.*//') || true

  # Mark completion (Python parses the actual log for all fields)
  ( flock -x 200
    echo "${graph},${r},${s},${algo},${status}" >> "$OUTCSV"
  ) 200>"$LOCKFILE"

  printf "  %-12s %-18s r=%s s=%-2s total=%sms peel=%sms status=%s\n" \
    "$algo" "$graph" "$r" "$s" "${total:-NA}" "${peel:-NA}" "$status"
  exit 0
fi

# ============ Step 0: Build ============
echo "============================================================="
echo "  V3 Phase Benchmark"
echo "  $(date)"
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
  echo "graph,r,s,algo,status" > "$OUTCSV"
fi
mkdir -p "$LOGDIR"

# ============ Generate jobs ============
GRAPHS=(com-dblp web-Stanford web-it-2004 dblp-core30 email-Eu-core com-youtube)
ALGOS="ST V3 V3_NP"

JOBFILE=$(mktemp /tmp/bench_v3_jobs.XXXXXX)

for graph in "${GRAPHS[@]}"; do
  if [ ! -f "graphs/${graph}.edges" ]; then
    echo "SKIP: graphs/${graph}.edges not found"
    continue
  fi

  for rs in "3 4" "3 5" "3 6" "3 8" "3 10" "4 5" "4 6" "5 6" "6 7" "6 8"; do
    r=${rs% *}; s=${rs#* }

    for algo in $ALGOS; do
      case $algo in
        ST)    env_var="PIVOTER_RUN_ST" ;;
        V3)    env_var="PIVOTER_RUN_REGION_V3" ;;
        V3_NP) env_var="PIVOTER_RUN_REGION_V3" ;;
      esac
      echo "--run $graph $r $s $algo $env_var" >> "$JOBFILE"
    done
  done
done

TOTAL=$(wc -l < "$JOBFILE")
EXISTING=$(tail -n +2 "$OUTCSV" 2>/dev/null | wc -l | tr -d ' ')
echo "Total jobs: $TOTAL, existing: $EXISTING, max parallel: $MAX_NPROC"
echo "Memory: launch<${MEM_LIMIT_GB}GB, kill>${MEM_KILL_GB}GB"
echo "Timeout per job: ${TIMEOUT}s"
echo "Output: $OUTCSV"
echo ""

# ============ Memory-aware job scheduler ============
# Strategy:
#   - Launch jobs as fast as possible (3s settle between launches)
#   - Don't launch if: jobs >= 32 OR memory usage >= 300GB
#   - If memory > 450GB: kill newest job
#     - If it was the ONLY job → it's genuinely OOM → mark OOM in CSV
#     - If there were multiple jobs → not its fault → re-queue for later
#   - This lets small jobs blast through 32-parallel, large jobs auto-throttle

get_used_mem_gb() {
  awk '/MemTotal/{t=$2} /MemAvailable/{a=$2} END{printf "%.0f", (t-a)/1024/1024}' /proc/meminfo 2>/dev/null || echo "0"
}

refresh_jobs() {
  local_pids=(); local_args=()
  for i in "${!JOB_PIDS[@]}"; do
    if kill -0 "${JOB_PIDS[$i]}" 2>/dev/null; then
      local_pids+=("${JOB_PIDS[$i]}")
      local_args+=("${JOB_ARGS[$i]}")
    fi
  done
  JOB_PIDS=("${local_pids[@]}")
  JOB_ARGS=("${local_args[@]}")
}

kill_newest_job() {
  local last=$(( ${#JOB_PIDS[@]} - 1 ))
  KILLED_ARGS="${JOB_ARGS[$last]}"
  local pid=${JOB_PIDS[$last]}
  # Kill bash wrapper — trap inside propagates to degeneracy_cliques child
  kill "$pid" 2>/dev/null
  # Wait for it to actually die (up to 10s)
  for _w in $(seq 10); do kill -0 "$pid" 2>/dev/null || break; sleep 1; done
  # Force kill if still alive
  kill -9 "$pid" 2>/dev/null
  unset 'JOB_PIDS[$last]'; unset 'JOB_ARGS[$last]'
  JOB_PIDS=("${JOB_PIDS[@]}"); JOB_ARGS=("${JOB_ARGS[@]}")
  echo "  [KILL] pid=$pid used=$(get_used_mem_gb)GB > ${MEM_KILL_GB}GB $(date +%H:%M:%S)"
}

write_oom() {
  local a="$1" g r s algo
  g=$(echo "$a" | awk '{print $2}'); r=$(echo "$a" | awk '{print $3}')
  s=$(echo "$a" | awk '{print $4}'); algo=$(echo "$a" | awk '{print $5}')
  ( flock -x 200
    echo "${g},${r},${s},${algo},OOM" >> "$OUTCSV"
  ) 200>"$LOCKFILE"
  echo "  [OOM] ${algo} ${g} r=${r} s=${s}"
}

# OOM check: kill newest jobs while memory > MEM_KILL_GB
oom_check() {
  while true; do
    local used=$(get_used_mem_gb)
    [ "$used" -lt "$MEM_KILL_GB" ] && return
    refresh_jobs
    [ ${#JOB_PIDS[@]} -eq 0 ] && return  # external memory pressure, nothing to kill
    local was_alone=$(( ${#JOB_PIDS[@]} == 1 ? 1 : 0 ))
    kill_newest_job
    if [ "$was_alone" -eq 1 ]; then
      # Only job running → genuinely OOM
      write_oom "$KILLED_ARGS"
    else
      # Multiple jobs → not its fault → re-queue
      RETRY_QUEUE+=("$KILLED_ARGS")
      echo "  → re-queued (was not alone, likely not its fault)"
    fi
    sleep 3
  done
}

SCRIPT="$(realpath "$0")"
LAUNCHED=0
declare -a JOB_PIDS=()
declare -a JOB_ARGS=()
declare -a RETRY_QUEUE=()

exec 3< "$JOBFILE"
FILE_DONE=false

while true; do
  # Pick next job: retry queue first, then job file
  jobargs=""
  if [ ${#RETRY_QUEUE[@]} -gt 0 ]; then
    jobargs="${RETRY_QUEUE[0]}"
    RETRY_QUEUE=("${RETRY_QUEUE[@]:1}")
  elif ! $FILE_DONE && IFS= read -r jobargs <&3; then
    :
  else
    FILE_DONE=true
    refresh_jobs
    [ ${#JOB_PIDS[@]} -eq 0 ] && [ ${#RETRY_QUEUE[@]} -eq 0 ] && break
    oom_check
    sleep $MEM_CHECK_INTERVAL
    continue
  fi

  # OOM protection before launching
  oom_check

  # Wait for: jobs < MAX_NPROC AND memory < MEM_LIMIT_GB
  while true; do
    refresh_jobs
    njobs=${#JOB_PIDS[@]}
    used=$(get_used_mem_gb)
    [ "$njobs" -lt "$MAX_NPROC" ] && [ "$used" -lt "$MEM_LIMIT_GB" ] && break
    [ "$used" -ge "$MEM_KILL_GB" ] && { oom_check; continue; }
    echo "  [wait] jobs=$njobs used=${used}GB $(date +%H:%M:%S)"
    sleep $MEM_CHECK_INTERVAL
  done

  # Launch (trap inside --run propagates kill to child)
  bash "$SCRIPT" $jobargs &
  JOB_PIDS+=($!)
  JOB_ARGS+=("$jobargs")
  LAUNCHED=$((LAUNCHED + 1))
  sleep $LAUNCH_SETTLE
done

exec 3<&-
rm -f "$JOBFILE" "$LOCKFILE"
echo ""
echo "Launched: $LAUNCHED jobs"

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV ($(wc -l < "$OUTCSV") rows)"
echo "  Logs: $LOGDIR/"
echo "  $(date)"
echo "============================================================="

# ============ Parse logs & generate summary ============
echo ""
echo "Parsing logs with Python..."
python3 bench_parse_logs.py
