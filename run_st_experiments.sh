#!/bin/bash
# Single-thread experiment: 4-way comparison
#   1. Reference (PIVOTER_RUN_REF=1, OMP_NUM_THREADS=1)
#   2. Optimized ST (PIVOTER_RUN_ST=1, OMP_NUM_THREADS=1)
#   3. Optimized T=1 (OMP_NUM_THREADS=1, existing parallel code at 1 thread)
#   4. Optimized T=16 (OMP_NUM_THREADS=16)
#
# Usage: ./run_st_experiments.sh [BINARY] [GRAPH_DIR] [RUNS]
#   BINARY   = path to degeneracy_cliques  (default: build/bin/degeneracy_cliques)
#   GRAPH_DIR = directory containing .edges files (default: graphs/)
#   RUNS     = number of runs per config (default: 3, reports median)

set -euo pipefail

BIN="${1:-build/bin/degeneracy_cliques}"
GRAPH_DIR="${2:-graphs/}"
RUNS="${3:-3}"

# Graphs to test (basenames without extension)
GRAPHS=(
    com-dblp
    com-youtube
    web-Google
    as-skitter
    soc-pokec
    wiki-Talk
)

# (r, s) configs with r in {1, 2}
declare -a CONFIGS=(
    "1 3"
    "1 4"
    "1 5"
    "2 3"
    "2 4"
    "2 5"
)

median() {
    # Read numbers from args, print median
    local -a sorted=($(printf '%s\n' "$@" | sort -n))
    local n=${#sorted[@]}
    echo "${sorted[$((n / 2))]}"
}

extract_time_ms() {
    # Extract "time: NNN ms" from output, return NNN (portable, no grep -oP)
    echo "$1" | sed -n 's/.*time: \([0-9]*\) ms.*/\1/p' | tail -1
}

echo "========================================"
echo "Single-Thread Experiment"
echo "Binary: $BIN"
echo "Graphs: ${GRAPHS[*]}"
echo "Runs per config: $RUNS"
echo "========================================"
printf "%-20s %-6s %-6s %12s %12s %12s %12s\n" \
    "Graph" "r" "s" "Ref(ms)" "ST(ms)" "T=1(ms)" "T=16(ms)"
echo "------------------------------------------------------------------------"

for graph in "${GRAPHS[@]}"; do
    graph_file="${GRAPH_DIR}/${graph}.edges"
    if [[ ! -f "$graph_file" ]]; then
        echo "SKIP: $graph_file not found"
        continue
    fi

    for config in "${CONFIGS[@]}"; do
        read -r r s <<< "$config"

        ref_times=()
        st_times=()
        t1_times=()
        t16_times=()

        for ((run = 0; run < RUNS; run++)); do
            # 1. Reference
            out=$(OMP_NUM_THREADS=1 PIVOTER_RUN_REF=1 "$BIN" "$graph_file" "$r" "$s" 2>&1 || true)
            t=$(extract_time_ms "$out")
            [[ -n "$t" ]] && ref_times+=("$t")

            # 2. Optimized ST
            out=$(OMP_NUM_THREADS=1 PIVOTER_RUN_ST=1 "$BIN" "$graph_file" "$r" "$s" 2>&1 || true)
            t=$(extract_time_ms "$out")
            [[ -n "$t" ]] && st_times+=("$t")

            # 3. Optimized T=1 (parallel code at 1 thread)
            out=$(OMP_NUM_THREADS=1 "$BIN" "$graph_file" "$r" "$s" 2>&1 || true)
            t=$(extract_time_ms "$out")
            [[ -n "$t" ]] && t1_times+=("$t")

            # 4. Optimized T=16
            out=$(OMP_NUM_THREADS=16 "$BIN" "$graph_file" "$r" "$s" 2>&1 || true)
            t=$(extract_time_ms "$out")
            [[ -n "$t" ]] && t16_times+=("$t")
        done

        ref_med=$(median "${ref_times[@]:-N/A}")
        st_med=$(median "${st_times[@]:-N/A}")
        t1_med=$(median "${t1_times[@]:-N/A}")
        t16_med=$(median "${t16_times[@]:-N/A}")

        printf "%-20s %-6s %-6s %12s %12s %12s %12s\n" \
            "$graph" "$r" "$s" "$ref_med" "$st_med" "$t1_med" "$t16_med"
    done
done

echo ""
echo "Interpretation:"
echo "  Ref → ST:   pure algorithmic speedup (no parallelism overhead)"
echo "  ST → T=1:   overhead from parallel infrastructure at 1 thread"
echo "  T=1 → T=16: parallel scaling factor"
