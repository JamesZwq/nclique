#!/usr/bin/env bash
# Launch all paper-r1 experiments in the background. Run on each server
# (tods, tods2) after `git pull`. The scripts auto-detect server name
# from $1 or hostname.
#
# Usage: ./run_all_experiments.sh [tods1|tods2]
#
# Each script picks its per-server graph set, builds the binary if
# needed, and writes results to paper_data/bench_*.csv. Per-job stdout
# goes to bench_*_logs/.
#
# Sequencing rationale:
#   * bench_par_sdct uses all OMP threads per job → must run alone.
#   * Other scripts use OMP_NUM_THREADS=1 + outer parallelism (many jobs
#     concurrently) → can run together as long as the shared memory
#     gates in bench_lib.py keep us under the per-server cap.
#
# Practical sequencing: run par-SDCT FIRST (sequential, finishes
# quickly), then launch the parallel-outer scripts in parallel.

set -u
SERVER="${1:-$(hostname)}"
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
mkdir -p paper_data run_logs

stamp() { date +%F_%H%M%S; }

echo "=== running on $SERVER ==="
echo

# Phase 1: parallel SDCT scaling — sequential outer, finishes first.
LOG=run_logs/par_sdct_$(stamp).log
echo "[1/6] bench_par_sdct (sequential outer) -> $LOG"
nohup python3 -u bench_par_sdct.py "$SERVER" > "$LOG" 2>&1 &
PID_PAR=$!
echo "    pid=$PID_PAR"
wait "$PID_PAR" || true

# Phase 2: parallel-outer scripts — launch all in parallel.
declare -a JOBS=(
    "bench_local_iterative.py:run_logs/local_$(stamp).log"
    "bench_ablation.py:run_logs/ablation_$(stamp).log"
    "bench_sort_options.py:run_logs/sort_$(stamp).log"
    "bench_hierarchy.py:run_logs/hierarchy_$(stamp).log"
    "bench_baselines.py:run_logs/baselines_$(stamp).log"
)

PIDS=()
for spec in "${JOBS[@]}"; do
    script=${spec%%:*}
    log=${spec##*:}
    echo "[ ] $script -> $log"
    nohup python3 -u "$script" "$SERVER" > "$log" 2>&1 &
    PIDS+=("$!")
done
echo
echo "Phase-2 PIDs: ${PIDS[*]}"
echo "Tail any: tail -f run_logs/*.log"
echo "Kill all: kill ${PIDS[*]}"

# Wait for all (so the parent stays alive under nohup).
for pid in "${PIDS[@]}"; do
    wait "$pid" || echo "  pid=$pid exited non-zero"
done
echo
echo "=== ALL DONE on $SERVER ==="
