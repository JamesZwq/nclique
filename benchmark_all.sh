#!/bin/bash
#
# Comprehensive benchmark script for all nucleus decomposition versions.
# Usage: bash benchmark_all.sh [output_dir]
#
# Runs all versions on multiple graphs with different r,s parameters.
# Parallel versions are tested with thread counts: 1, 2, 4, 8, 16, 32, 64.
# Records wall-clock time and peak RSS memory.
#

set -euo pipefail

# Build first
echo "Building..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build --target degeneracy_cliques -j $(nproc) 2>&1 | tail -3
echo "Build complete."
echo ""

BIN="./build/bin/degeneracy_cliques"
OUTDIR="${1:-benchmark_results}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="${OUTDIR}/${TIMESTAMP}"
mkdir -p "$OUTDIR"

DATASETS=(
    "/data/wenqianz/com-dblp.edges"
    "/data/wenqianz/com-youtube.edges"
    "/data/wenqianz/web-Google.edges"
    "/data/wenqianz/soc-pokec-relationships.edges"
)

THREAD_COUNTS=(1 2 4 8 16 32 64)

# =========================================================================
# Helper: run one experiment, capture time + peak RSS
# Args: $1=label $2=env_vars $3=threads $4=graph $5=r $6=s $7=outfile
# =========================================================================
run_one() {
    local label="$1"
    local env_vars="$2"
    local threads="$3"
    local graph="$4"
    local r_val="$5"
    local s_val="$6"
    local outfile="$7"
    local gname
    gname=$(basename "$graph" .edges)

    echo -n "  Running: ${label} threads=${threads} ${gname} r=${r_val} s=${s_val} ... "

    # Use /usr/bin/time for peak RSS (GNU time -v or BSD time -l)
    local tmplog
    tmplog=$(mktemp /tmp/pivoter_bench.XXXXXX)

    # Set timeout to 30 minutes
    local cmd="OMP_NUM_THREADS=${threads} ${env_vars} timeout 1800 ${BIN} ${graph} ${r_val} ${s_val} degen"

    # Capture wall time and peak RSS via /usr/bin/time
    { /usr/bin/time -v bash -c "$cmd" ; } > "$tmplog" 2>&1 || true

    # Extract peeling time (the "time: XXX ms" line from the algorithm)
    local peel_time
    peel_time=$(grep -oP 'time: \K[0-9]+' "$tmplog" | tail -1 || echo "N/A")

    # Extract peak RSS from GNU time output
    local peak_rss_kb
    peak_rss_kb=$(grep -oP 'Maximum resident set size.*?: \K[0-9]+' "$tmplog" || echo "N/A")

    # Extract total wall time from the "took:" line
    local total_time
    total_time=$(grep -oP '(?:ST r|Reference r|Local H-index|NucleusCoreDecomposition) .* took: \K[0-9.]+' "$tmplog" | tail -1 || echo "N/A")

    # Extract Cases info if available
    local cases_info
    cases_info=$(grep -oP 'Cases: .*' "$tmplog" | tail -1 || echo "")

    # Write CSV row
    echo "${gname},${r_val},${s_val},${label},${threads},${peel_time},${total_time},${peak_rss_kb},${cases_info}" >> "$outfile"

    echo "peel=${peel_time}ms rss=${peak_rss_kb}kB"
    rm -f "$tmplog"
}

# =========================================================================
# Main benchmark loop
# =========================================================================

echo "=============================================="
echo " Pivoter Benchmark Suite"
echo " Output: ${OUTDIR}"
echo " Date: $(date)"
echo "=============================================="
echo ""

# CSV header
CSVFILE="${OUTDIR}/results.csv"
echo "graph,r,s,version,threads,peel_time_ms,total_time_ms,peak_rss_kb,cases_info" > "$CSVFILE"

# ===================== r=1 benchmarks =====================
echo "====== r=1 benchmarks ======"

for graph in "${DATASETS[@]}"; do
    for s_val in 3 4 5; do
        gname=$(basename "$graph" .edges)
        echo "--- ${gname} r=1 s=${s_val} ---"

        # Correct (single-thread)
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 1 "$s_val" "$CSVFILE"

        # ST peeling (single-thread)
        run_one "ST" "PIVOTER_RUN_ST=1" 1 "$graph" 1 "$s_val" "$CSVFILE"

        # Local V1 (single-thread)
        run_one "Local_V1" "PIVOTER_RUN_LOCAL=1" 1 "$graph" 1 "$s_val" "$CSVFILE"

        # Local V2 (single-thread)
        run_one "Local_V2" "PIVOTER_RUN_LOCAL_V2=1" 1 "$graph" 1 "$s_val" "$CSVFILE"

        # Local V3 parallel — sweep threads
        for t in "${THREAD_COUNTS[@]}"; do
            run_one "Local_V3" "PIVOTER_RUN_LOCAL_V3=1" "$t" "$graph" 1 "$s_val" "$CSVFILE"
        done

        # Local V4 parallel — sweep threads
        for t in "${THREAD_COUNTS[@]}"; do
            run_one "Local_V4" "PIVOTER_RUN_LOCAL_V4=1" "$t" "$graph" 1 "$s_val" "$CSVFILE"
        done

        echo ""
    done
done

# ===================== r=2 benchmarks =====================
echo "====== r=2 benchmarks ======"

for graph in "${DATASETS[@]}"; do
    for s_val in 3 4 5; do
        gname=$(basename "$graph" .edges)
        echo "--- ${gname} r=2 s=${s_val} ---"

        # Correct (single-thread)
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 2 "$s_val" "$CSVFILE"

        # ST peeling (single-thread)
        run_one "ST" "PIVOTER_RUN_ST=1" 1 "$graph" 2 "$s_val" "$CSVFILE"

        # Plus parallel (default r=2) — sweep threads
        for t in "${THREAD_COUNTS[@]}"; do
            run_one "Plus_Parallel" "" "$t" "$graph" 2 "$s_val" "$CSVFILE"
        done

        echo ""
    done
done

# ===================== r=3 benchmarks =====================
echo "====== r=3 benchmarks ======"

for graph in "${DATASETS[@]}"; do
    for s_val in 4 5; do
        gname=$(basename "$graph" .edges)
        echo "--- ${gname} r=3 s=${s_val} ---"

        # Correct (single-thread)
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 3 "$s_val" "$CSVFILE"

        # ST (single-thread)
        run_one "ST" "PIVOTER_RUN_ST=1" 1 "$graph" 3 "$s_val" "$CSVFILE"

        # Parallel (default r>=3) — sweep threads
        for t in "${THREAD_COUNTS[@]}"; do
            run_one "RClique_Parallel" "" "$t" "$graph" 3 "$s_val" "$CSVFILE"
        done

        echo ""
    done
done

# ===================== r=4 benchmarks =====================
echo "====== r=4 benchmarks ======"

for graph in "${DATASETS[@]}"; do
    for s_val in 5; do
        gname=$(basename "$graph" .edges)
        echo "--- ${gname} r=4 s=${s_val} ---"

        # Correct (single-thread)
        run_one "Correct" "PIVOTER_RUN_REF=1" 1 "$graph" 4 "$s_val" "$CSVFILE"

        # ST (single-thread)
        run_one "ST" "PIVOTER_RUN_ST=1" 1 "$graph" 4 "$s_val" "$CSVFILE"

        # Parallel (default r>=3) — sweep threads
        for t in "${THREAD_COUNTS[@]}"; do
            run_one "RClique_Parallel" "" "$t" "$graph" 4 "$s_val" "$CSVFILE"
        done

        echo ""
    done
done

echo "=============================================="
echo " Benchmark complete!"
echo " Results: ${CSVFILE}"
echo "=============================================="
