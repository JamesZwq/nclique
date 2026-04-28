#!/usr/bin/env bash
# Relaunch only the experiments affected by the V3 r=1 dispatch bug:
#   * bench_ablation
#   * bench_local_iterative
#   * bench_sort_options
# All three need the new degeneracy_cliques binary (built from commit
# 88e047c or later). par_sdct, hierarchy, baselines are NOT relaunched.

set -u
SERVER="${1:-$(hostname)}"
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
mkdir -p run_logs
stamp() { date +%F_%H%M%S; }

echo "=== relaunch on $SERVER ==="

declare -a JOBS=(
    "bench_local_iterative.py:run_logs/local_$(stamp).log"
    "bench_ablation.py:run_logs/ablation_$(stamp).log"
    "bench_sort_options.py:run_logs/sort_$(stamp).log"
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
echo "PIDs: ${PIDS[*]}"
echo "Tail any: tail -f run_logs/*.log"
echo "Kill all: kill ${PIDS[*]}"

for pid in "${PIDS[@]}"; do
    wait "$pid" || echo "  pid=$pid exited non-zero"
done
echo
echo "=== ALL DONE on $SERVER ==="
