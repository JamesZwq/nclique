#!/usr/bin/env bash
# Sequential relaunch: run each script alone (one at a time) to avoid the
# 3-way memory contention that caused 30-50% OOM_KILL_NEWEST in the
# parallel-outer launcher. Each script's harness already does outer
# parallelism within its budget; running them serially just removes
# inter-script competition.
#
# The harness's resume support (load_done) only treats status=OK as
# done, so any OOM_KILL_NEWEST / SIGSEGV / TIMEOUT / PARSE_FAIL cells
# in the existing CSVs will be re-attempted automatically.
#
# Order picked to fail fast on issues:
#   1. par_sdct  — recently fixed (leaves overflow → double); fastest
#                  to verify the fix works.
#   2. sort_options — most missing data (only 53/65 OK on tods1).
#   3. hierarchy   — only 8 pairs collected; rerun to expand coverage.
#   4. ablation    — fill OOM_KILL_NEWEST cells (V2 SIGSEGVs will retry
#                    and fail again — that's expected, ~20 min).
#   5. local_iter  — Phase A OOM cells + Phase B re-attempt for any OOM.
#                    LONG (hours) — last so the others land first.
#
# Usage: ./relaunch_sequential.sh [tods1|tods2]

set -u
SERVER="${1:-$(hostname)}"
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
mkdir -p run_logs

stamp() { date +%F_%H%M%S; }

declare -a SCRIPTS=(
    "bench_par_sdct.py"
    "bench_sort_options.py"
    "bench_hierarchy.py"
    "bench_ablation.py"
    "bench_local_iterative.py"
)

for script in "${SCRIPTS[@]}"; do
    log="run_logs/${script%.py}_$(stamp).log"
    echo "=================================================="
    echo "[$(date '+%F %T')] $script -> $log"
    echo "=================================================="
    python3 -u "$script" "$SERVER" 2>&1 | tee "$log"
    echo "[$(date '+%F %T')] $script DONE"
    echo
done

echo "=== ALL SEQUENTIAL DONE on $SERVER ==="
