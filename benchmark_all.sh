#!/bin/bash
#
# Comprehensive benchmark script for all nucleus decomposition versions.
# Usage: bash benchmark_all.sh [--verify] [output_dir]
#
#   --verify   First run correctness verification (PIVOTER_COMPARE=1),
#              then run performance benchmarks.
#   (default)  Skip verification, run performance benchmarks only.
#
# All program output (memory, timing at each phase) is saved per-run in logs/.
# Single-thread only. Records wall-clock time and peak RSS memory.
#

set -euo pipefail

# =========================================================================
# Parse arguments
# =========================================================================
VERIFY=0
POSITIONAL=()
for arg in "$@"; do
    case "$arg" in
        --verify) VERIFY=1 ;;
        *) POSITIONAL+=("$arg") ;;
    esac
done

# Project root on the server
PROJDIR="${PROJDIR:-/home/wenqianz/nclique}"
cd "$PROJDIR"

# Pull latest and build
echo "Pulling latest code..."
git pull 2>&1 | tail -3
echo ""
echo "Building..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build --target degeneracy_cliques -j 12 2>&1 | tail -3
echo "Build complete."
echo ""

BIN="./build/bin/degeneracy_cliques"
OUTDIR="${POSITIONAL[0]:-benchmark_results}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="${OUTDIR}/${TIMESTAMP}"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$LOGDIR"

DATASETS=(
    "/data/wenqianz/com-dblp.edges"
    "/data/wenqianz/com-youtube.edges"
    "/data/wenqianz/web-Google.edges"
    "/data/wenqianz/soc-pokec-relationships.edges"
)

# =========================================================================
# Helper: run one experiment, capture time + peak RSS + full program output
# Args: $1=label $2=env_vars $3=threads $4=graph $5=r $6=s $7=csvfile
# =========================================================================
run_one() {
    local label="$1"
    local env_vars="$2"
    local threads="$3"
    local graph="$4"
    local r_val="$5"
    local s_val="$6"
    local csvfile="$7"
    local gname
    gname=$(basename "$graph" .edges)

    echo -n "  Running: ${label} ${gname} r=${r_val} s=${s_val} ... "

    # Per-run log file: stores ALL program output
    local logfile="${LOGDIR}/${gname}_r${r_val}_s${s_val}_${label}.log"

    local cmd="OMP_NUM_THREADS=${threads} ${env_vars} timeout 10800 ${BIN} ${graph} ${r_val} ${s_val} degen"

    # Capture everything: program stdout/stderr + /usr/bin/time stats
    { /usr/bin/time -v bash -c "$cmd" ; } > "$logfile" 2>&1 || true

    # Extract peeling time (the "time: XXX ms" line from the algorithm)
    local peel_time
    peel_time=$(grep -oP 'time: \K[0-9]+' "$logfile" | tail -1 || echo "N/A")

    # Extract peak RSS from GNU time output
    local peak_rss_kb
    peak_rss_kb=$(grep -oP 'Maximum resident set size.*?: \K[0-9]+' "$logfile" || echo "N/A")

    # Extract total wall time from the "took:" line
    local total_time
    total_time=$(grep -oP '(?:ST r|ST_V2 r|Reference r|Local H-index|NucleusCoreDecomposition) .* took: \K[0-9.]+' "$logfile" | tail -1 || echo "N/A")

    # Extract Cases info if available
    local cases_info
    cases_info=$(grep -oP 'Cases: .*' "$logfile" | tail -1 || echo "")

    # Write CSV row
    echo "${gname},${r_val},${s_val},${label},${threads},${peel_time},${total_time},${peak_rss_kb},${cases_info}" >> "$csvfile"

    echo "peel=${peel_time}ms rss=${peak_rss_kb}kB"
}

# =========================================================================
# Helper: run correctness verification for one (graph, r, s, version)
# Args: $1=label $2=env_vars $3=graph $4=r $5=s
# Returns 0 if correct, 1 if mismatch
# =========================================================================
verify_one() {
    local label="$1"
    local env_vars="$2"
    local graph="$3"
    local r_val="$4"
    local s_val="$5"
    local gname
    gname=$(basename "$graph" .edges)

    echo -n "  Verify: ${label} ${gname} r=${r_val} s=${s_val} ... "

    local logfile="${LOGDIR}/${gname}_r${r_val}_s${s_val}_${label}_verify.log"
    local cmd="OMP_NUM_THREADS=1 PIVOTER_COMPARE=1 ${env_vars} timeout 10800 ${BIN} ${graph} ${r_val} ${s_val} degen"

    { /usr/bin/time -v bash -c "$cmd" ; } > "$logfile" 2>&1 || true

    if grep -q '✓.*correctness verified' "$logfile"; then
        echo "✓ PASS"
        return 0
    elif grep -q '✗.*MISMATCH' "$logfile"; then
        echo "✗ FAIL (see ${logfile})"
        return 1
    else
        echo "? UNKNOWN (see ${logfile})"
        return 1
    fi
}

# =========================================================================
# Algorithm tables: (env_var, label) pairs per r-value
# =========================================================================

R1_VARIANTS=(
    "PIVOTER_RUN_ST=1        ST"
    "PIVOTER_RUN_ST_V2=1     ST_V2"
    "PIVOTER_RUN_ONDEMAND=1  OnDemand"
    "PIVOTER_RUN_INTERLEAVED=1 Interleaved"
    "PIVOTER_RUN_INTERLEAVED_V2=1 InterleavedV2"
    "PIVOTER_RUN_LOCAL=1     Local_V1"
    "PIVOTER_RUN_LOCAL_V2=1  Local_V2"
    "PIVOTER_RUN_LOCAL_V3=1  Local_V3"
    "PIVOTER_RUN_LOCAL_V4=1  Local_V4"
)

R2_VARIANTS=(
    "PIVOTER_RUN_ST=1        ST"
    "PIVOTER_RUN_ST_V4=1     ST_V4"
    "PIVOTER_RUN_R2_ONDEMAND=1 R2_OnDemand"
    "PIVOTER_RUN_R2_TREEFREE=1 R2_TreeFree"
    "PIVOTER_RUN_R2_TREEFREE_V2=1 R2_TreeFreeV2"
)

R3_VARIANTS=(
    "PIVOTER_RUN_ST=1        ST"
    "PIVOTER_RUN_ST_V4=1     ST_V4"
    "PIVOTER_RUN_ST_V10=1    ST_V10"
    "PIVOTER_RUN_ST_V11=1    ST_V11"
    "PIVOTER_RUN_ST_V12=1    ST_V12"
    "PIVOTER_RUN_ST_V13=1    ST_V13"
)

# =========================================================================
# Main
# =========================================================================

echo "=============================================="
echo " Pivoter Benchmark Suite"
echo " Output : ${OUTDIR}"
echo " Verify : $([ $VERIFY -eq 1 ] && echo YES || echo NO)"
echo " Date   : $(date)"
echo "=============================================="
echo ""

CSVFILE="${OUTDIR}/results.csv"
echo "graph,r,s,version,threads,peel_time_ms,total_time_ms,peak_rss_kb,cases_info" > "$CSVFILE"

VERIFY_FAIL=0

# =========================================================================
# Phase A: Correctness verification (optional)
# =========================================================================
if [ $VERIFY -eq 1 ]; then
    echo "========== Phase A: Correctness Verification =========="
    echo ""

    for graph in "${DATASETS[@]}"; do
        gname=$(basename "$graph" .edges)

        # r=1
        for s_val in 3 4 5; do
            echo "--- ${gname} r=1 s=${s_val} ---"
            for entry in "${R1_VARIANTS[@]}"; do
                env_var=$(echo "$entry" | awk '{print $1}')
                label=$(echo "$entry" | awk '{print $2}')
                verify_one "$label" "$env_var" "$graph" 1 "$s_val" || VERIFY_FAIL=1
            done
            echo ""
        done

        # r=2
        for s_val in 3 4 5; do
            echo "--- ${gname} r=2 s=${s_val} ---"
            for entry in "${R2_VARIANTS[@]}"; do
                env_var=$(echo "$entry" | awk '{print $1}')
                label=$(echo "$entry" | awk '{print $2}')
                verify_one "$label" "$env_var" "$graph" 2 "$s_val" || VERIFY_FAIL=1
            done
            echo ""
        done

        # r=3
        for s_val in 4 5; do
            echo "--- ${gname} r=3 s=${s_val} ---"
            for entry in "${R3_VARIANTS[@]}"; do
                env_var=$(echo "$entry" | awk '{print $1}')
                label=$(echo "$entry" | awk '{print $2}')
                verify_one "$label" "$env_var" "$graph" 3 "$s_val" || VERIFY_FAIL=1
            done
            echo ""
        done

        # r=4
        for s_val in 5; do
            echo "--- ${gname} r=4 s=${s_val} ---"
            for entry in "${R3_VARIANTS[@]}"; do
                env_var=$(echo "$entry" | awk '{print $1}')
                label=$(echo "$entry" | awk '{print $2}')
                verify_one "$label" "$env_var" "$graph" 4 "$s_val" || VERIFY_FAIL=1
            done
            echo ""
        done
    done

    if [ $VERIFY_FAIL -ne 0 ]; then
        echo "!!! CORRECTNESS FAILURES DETECTED — check logs in ${LOGDIR} !!!"
        echo "Aborting benchmark."
        exit 1
    fi
    echo "========== All correctness checks passed =========="
    echo ""
fi

# =========================================================================
# Phase B: Performance benchmarks
# =========================================================================
echo "========== Phase B: Performance Benchmarks =========="
echo ""

for graph in "${DATASETS[@]}"; do
    gname=$(basename "$graph" .edges)

    # r=1
    for s_val in 3 4 5; do
        echo "--- ${gname} r=1 s=${s_val} ---"
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 1 "$s_val" "$CSVFILE"
        for entry in "${R1_VARIANTS[@]}"; do
            env_var=$(echo "$entry" | awk '{print $1}')
            label=$(echo "$entry" | awk '{print $2}')
            run_one "$label" "$env_var" 1 "$graph" 1 "$s_val" "$CSVFILE"
        done
        echo ""
    done

    # r=2
    for s_val in 3 4 5; do
        echo "--- ${gname} r=2 s=${s_val} ---"
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 2 "$s_val" "$CSVFILE"
        for entry in "${R2_VARIANTS[@]}"; do
            env_var=$(echo "$entry" | awk '{print $1}')
            label=$(echo "$entry" | awk '{print $2}')
            run_one "$label" "$env_var" 1 "$graph" 2 "$s_val" "$CSVFILE"
        done
        echo ""
    done

    # r=3
    for s_val in 4 5; do
        echo "--- ${gname} r=3 s=${s_val} ---"
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 3 "$s_val" "$CSVFILE"
        for entry in "${R3_VARIANTS[@]}"; do
            env_var=$(echo "$entry" | awk '{print $1}')
            label=$(echo "$entry" | awk '{print $2}')
            run_one "$label" "$env_var" 1 "$graph" 3 "$s_val" "$CSVFILE"
        done
        echo ""
    done

    # r=4
    for s_val in 5; do
        echo "--- ${gname} r=4 s=${s_val} ---"
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 4 "$s_val" "$CSVFILE"
        for entry in "${R3_VARIANTS[@]}"; do
            env_var=$(echo "$entry" | awk '{print $1}')
            label=$(echo "$entry" | awk '{print $2}')
            run_one "$label" "$env_var" 1 "$graph" 4 "$s_val" "$CSVFILE"
        done
        echo ""
    done
done

echo "=============================================="
echo " Benchmark complete!"
echo " Results CSV : ${CSVFILE}"
echo " Full logs   : ${LOGDIR}/"
echo "=============================================="
