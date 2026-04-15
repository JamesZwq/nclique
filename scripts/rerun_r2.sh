#!/bin/bash
# Rerun all R=2 experiments with CSR-release fix
# Wait for benchmark_all to finish, then rebuild and rerun

cd "$(dirname "$0")/.."

BIN="./build/bin/degeneracy_cliques"
TIMEOUT=3600
OUTCSV="benchmark_r2_rerun.csv"
LOGDIR="bench_logs_r2_rerun"
LOCKFILE="/tmp/bench_r2_rerun.lock"
NPROC=4

# Wait for benchmark_all to finish
echo "Waiting for benchmark_all to finish..."
while pgrep -f "benchmark_all" > /dev/null 2>&1; do
    sleep 60
done
echo "benchmark_all done."

# Rebuild
echo "Rebuilding..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build -j 12 --target degeneracy_cliques 2>&1 | tail -3
echo "Build done."

# Setup
mkdir -p "$LOGDIR"
if [ ! -f "$OUTCSV" ] || ! head -1 "$OUTCSV" | grep -q 'graph,r,s'; then
    echo "graph,r,s,algorithm,time_ms,memory_kB,status" > "$OUTCSV"
fi

GRAPHS=(com-youtube web-Stanford web-it-2004)
MAX_CLIQUE=(17 61 430)

JOBFILE=$(mktemp /tmp/bench_r2_jobs.XXXXXX)

for gi in "${!GRAPHS[@]}"; do
    graph="${GRAPHS[$gi]}"
    omega="${MAX_CLIQUE[$gi]}"
    [ ! -f "graphs/${graph}.edges" ] && continue

    for s in $(seq 3 1 "$omega"); do
        echo "$graph 2 $s Ours_DCLP PIVOTER_RUN_R2_DCLP" >> "$JOBFILE"
        echo "$graph 2 $s REF_R2 PIVOTER_NOOPT" >> "$JOBFILE"
    done
done

TOTAL=$(wc -l < "$JOBFILE" | tr -d ' ')
echo "Total R=2 jobs: $TOTAL"
echo "Running with $NPROC parallel threads..."

run_one() {
    local graph=$1 r=$2 s=$3 algo=$4 env_var=$5
    local logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo}.log"

    if grep -q "^${graph},${r},${s},${algo}," "$OUTCSV" 2>/dev/null; then
        return 0
    fi

    local result="" exit_code=0
    result=$(env "${env_var}=1" timeout "${TIMEOUT}s" "$BIN" "graphs/${graph}.edges" "$r" "$s" 2>&1) || exit_code=$?
    echo "$result" > "$logfile"

    local took="" status="OK"
    if [ "$exit_code" -eq 124 ]; then
        took="TIMEOUT"; status="TIMEOUT"
    else
        took=$(echo "$result" | grep -oP 'NucleusCoreDecomposition took: \K[0-9.]+(?= ms)' || true)
        [ -z "$took" ] && took=$(echo "$result" | grep "^time:" | tail -1 | awk '{print $2}' || true)
        [ -z "$took" ] && { took="ERROR"; status="ERROR(exit=${exit_code})"; }
    fi

    local mem
    mem=$(echo "$result" | grep -oP 'Final Memory: \s*\K[0-9]+' || true)
    [ -z "$mem" ] && mem="N/A"

    (
        flock -x 200
        echo "${graph},${r},${s},${algo},${took},${mem},${status}" >> "$OUTCSV"
    ) 200>"$LOCKFILE"

    printf "    %-12s r=%s s=%-3s %10sms (%s) mem=%s kB\n" "$algo" "$r" "$s" "$took" "$status" "$mem"
}

export -f run_one
export BIN TIMEOUT OUTCSV LOGDIR LOCKFILE

cat "$JOBFILE" | while read graph r s algo env_var; do
    echo "$graph $r $s $algo $env_var"
done | xargs -P "$NPROC" -L 1 bash -c 'run_one $@' _

rm -f "$JOBFILE" "$LOCKFILE"
echo ""
echo "Done. Results: $OUTCSV"
echo "Logs: $LOGDIR/"
