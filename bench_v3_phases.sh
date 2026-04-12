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
  total=$(echo "$result" | grep 'NucleusCoreDecomposition took:' | awk '{print $3}') || true

  # Extract memory
  graph_mem=$(echo "$result" | grep "Graph Memory" | awk -F: '{print $NF}' | tr -d ' kB') || true
  index_mem=$(echo "$result" | grep -E "Other Index Memory|Other index" | awk -F: '{print $NF}' | tr -d ' kB') || true
  final_mem=$(echo "$result" | grep "Final Memory" | awk -F: '{print $NF}' | tr -d ' kB') || true

  # V3-specific phases
  sdct_time=$(echo "$result" | grep 'SDCT_MaxClique took:' | awk '{print $3}') || true
  maxcliq_time=$(echo "$result" | grep 'MaxCliqEnum (V3) took:' | awk '{print $4}') || true
  merge_time=$(echo "$result" | grep 'r-Mergeable classification:' | awk '{print $3}') || true
  cpi_time=$(echo "$result" | grep "CPI counting time:" | awk -F: '{print $NF}' | tr -d ' ms') || true
  pathinfo_time=$(echo "$result" | grep "PathInfo build time:" | awk -F: '{print $NF}' | tr -d ' ms') || true
  peel_time=$(echo "$result" | grep "Peeling time:" | awk -F: '{print $NF}' | tr -d ' ms') || true
  v3_total=$(echo "$result" | grep "Total time:" | awk -F: '{print $NF}' | tr -d ' ms') || true

  # V3-specific stats
  fully_merge=$(echo "$result" | grep "Fully mergeable" | awk '{print $4}') || true
  remaining=$(echo "$result" | grep "Remaining regions:" | awk -F: '{print $NF}' | tr -d ' ') || true
  classes=$(echo "$result" | grep "Overlap classes:" | awk -F: '{print $NF}' | tr -d ' ') || true
  tuples=$(echo "$result" | grep "r-tuples:" | awk -F: '{print $NF}' | tr -d ' ') || true
  recursive=$(echo "$result" | grep "Total recursive calls:" | awk -F: '{print $NF}' | tr -d ' ') || true
  max_core=$(echo "$result" | grep "Max core:" | head -1 | awk -F: '{print $NF}' | tr -d ' ') || true
  rcliques=$(echo "$result" | grep "r-cliques:" | head -1 | awk -F: '{print $NF}' | tr -d ' ') || true

  # V2-specific
  v2_enum=$(echo "$result" | grep "Enumeration time:" | awk -F: '{print $NF}' | tr -d ' ms') || true
  v2_stuples=$(echo "$result" | grep "s-tuples:" | head -1 | awk -F: '{print $NF}' | tr -d ' ') || true
  v2_peel=$(echo "$result" | grep "Peeling time:" | awk -F: '{print $NF}' | tr -d ' ms') || true

  # ST-specific
  st_fused=$(echo "$result" | grep 'fused build+counting' | awk '{print $4}') || true
  st_bk=$(echo "$result" | grep "BK:" | awk -F: '{print $NF}' | tr -d ' ') || true

  (
    flock -x 200
    echo "${graph},${r},${s},${algo},${status},${total:-NA},${graph_mem:-NA},${index_mem:-NA},${final_mem:-NA},${sdct_time:-NA},${maxcliq_time:-NA},${merge_time:-NA},${cpi_time:-NA},${pathinfo_time:-NA},${peel_time:-NA},${v3_total:-NA},${fully_merge:-NA},${remaining:-NA},${classes:-NA},${tuples:-NA},${recursive:-NA},${max_core:-NA},${rcliques:-NA},${v2_enum:-NA},${v2_stuples:-NA},${v2_peel:-NA},${st_fused:-NA},${st_bk:-NA}" >> "$OUTCSV"
  ) 200>"$LOCKFILE"

  printf "  %-12s %-18s r=%s s=%-2s total=%sms peel=%sms status=%s\n" \
    "$algo" "$graph" "$r" "$s" "${total:-NA}" "${peel_time:-NA}" "$status"
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
  echo "graph,r,s,algo,status,total_ms,graph_mem_kB,index_mem_kB,final_mem_kB,sdct_ms,maxcliq_ms,merge_ms,cpi_ms,pathinfo_ms,peel_ms,v3_total_ms,fully_merge,remaining,classes,tuples,recursive_calls,max_core,rcliques,v2_enum_ms,v2_stuples,v2_peel_ms,st_fused_ms,st_bk_ms" > "$OUTCSV"
fi
mkdir -p "$LOGDIR"

# ============ Generate jobs ============
GRAPHS=(com-dblp web-Stanford web-it-2004 dblp-core30 email-Eu-core com-youtube)
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
  kill "$pid" 2>/dev/null; wait "$pid" 2>/dev/null
  unset 'JOB_PIDS[$last]'; unset 'JOB_ARGS[$last]'
  JOB_PIDS=("${JOB_PIDS[@]}"); JOB_ARGS=("${JOB_ARGS[@]}")
  echo "  [KILL] pid=$pid used=$(get_used_mem_gb)GB > ${MEM_KILL_GB}GB $(date +%H:%M:%S)"
}

write_oom() {
  local a="$1" g r s algo
  g=$(echo "$a" | awk '{print $2}'); r=$(echo "$a" | awk '{print $3}')
  s=$(echo "$a" | awk '{print $4}'); algo=$(echo "$a" | awk '{print $5}')
  ( flock -x 200
    echo "${g},${r},${s},${algo},OOM,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> "$OUTCSV"
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

  # Launch
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

# ============ Summary table ============
echo ""
echo "=== Summary: Total Time (ms) ==="
echo ""
printf "%-18s %3s %3s | %10s %10s %10s\n" "Graph" "r" "s" "ST" "V2" "V3"
echo "---------------------------------------------------------------------"

tail -n +2 "$OUTCSV" | sort -t',' -k1,1 -k2,2n -k3,3n | awk -F',' '
{
  key=$1","$2","$3
  if ($4=="ST") st[key]=$6
  if ($4=="V2") v2[key]=$6
  if ($4=="V3") v3[key]=$6
  seen[key]=1
}
END {
  for (key in seen) {
    split(key, a, ",")
    printf "%-18s %3s %3s | %10s %10s %10s\n",
      a[1], a[2], a[3],
      (st[key]?st[key]:"NA"), (v2[key]?v2[key]:"NA"), (v3[key]?v3[key]:"NA")
  }
}' | sort

echo ""
echo "=== Summary: Memory (MB) ==="
echo ""
printf "%-18s %3s %3s | %8s %8s %8s\n" "Graph" "r" "s" "ST" "V2" "V3"
echo "---------------------------------------------------------------------"

tail -n +2 "$OUTCSV" | sort -t',' -k1,1 -k2,2n -k3,3n | awk -F',' '
{
  key=$1","$2","$3
  if ($4=="ST") st[key]=$9
  if ($4=="V2") v2[key]=$9
  if ($4=="V3") v3[key]=$9
  seen[key]=1
}
END {
  for (key in seen) {
    split(key, a, ",")
    st_mb = (st[key] && st[key]!="NA") ? sprintf("%.0f", st[key]/1024) : "NA"
    v2_mb = (v2[key] && v2[key]!="NA") ? sprintf("%.0f", v2[key]/1024) : "NA"
    v3_mb = (v3[key] && v3[key]!="NA") ? sprintf("%.0f", v3[key]/1024) : "NA"
    printf "%-18s %3s %3s | %8s %8s %8s\n",
      a[1], a[2], a[3], st_mb, v2_mb, v3_mb
  }
}' | sort
