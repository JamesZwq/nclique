#!/bin/bash
# Full experiment script for tods2 server
# Usage: nohup bash run_full_experiments.sh > experiment_run.log 2>&1 &

PROJECT_DIR="$HOME/nclique_tmp"
BIN="$PROJECT_DIR/build/bin/degeneracy_cliques"
LOGDIR="$PROJECT_DIR/experiment_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR"

GRAPHS=(
    "/data/wenqianz/com-dblp.edges"
    "/data/wenqianz/web-Google.edges"
    "/data/wenqianz/web-Stanford.edges"
    "/data/wenqianz/com-youtube.edges"
    "/data/wenqianz/soc-pokec-relationships.edges"
)
GRAPH_NAMES=(
    "com-dblp"
    "web-Google"
    "web-Stanford"
    "com-youtube"
    "soc-pokec"
)

RS_PAIRS=("1 3" "1 4" "1 5" "2 3" "2 4" "2 5" "3 4" "3 5" "4 5")

THREADS=16
TIMEOUT_SEC=3600

# Helper: extract "took: NNN.NN ms" value after a keyword
extract_took() {
    echo "$1" | grep "$2" | grep "took:" | sed 's/.*took: \([0-9.]*\) ms.*/\1/' | head -1
}

echo "============================================================"
echo "  FULL EXPERIMENT SUITE"
echo "  Server: $(hostname), Cores: $(nproc), Threads: $THREADS"
echo "  Date: $(date)"
echo "  Log directory: $LOGDIR"
echo "============================================================"

# ============================================================
# Part 0: Build
# ============================================================
echo ""
echo ">>> Building project..."
cd "$PROJECT_DIR"
rm -rf build
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build -j$(nproc) --target degeneracy_cliques 2>&1 | tail -5
if [ ! -x "$BIN" ]; then
    echo "FATAL: Build failed, $BIN not found"
    exit 1
fi
echo "Build successful: $BIN"
echo ""

# ============================================================
# Part 1: Correctness verification
# ============================================================
VERIFY_LOG="$LOGDIR/correctness_verification.log"
echo "============================================================" | tee "$VERIFY_LOG"
echo "  PART 1: Correctness Verification" | tee -a "$VERIFY_LOG"
echo "============================================================" | tee -a "$VERIFY_LOG"

VERIFY_PASS=0
VERIFY_FAIL=0
VERIFY_SKIP=0

for gi in "${!GRAPHS[@]}"; do
    graph="${GRAPHS[$gi]}"
    gname="${GRAPH_NAMES[$gi]}"
    [ ! -f "$graph" ] && echo "SKIP: $graph not found" | tee -a "$VERIFY_LOG" && continue

    for rs in "${RS_PAIRS[@]}"; do
        r=$(echo $rs | cut -d' ' -f1)
        s=$(echo $rs | cut -d' ' -f2)
        echo -n "Verify $gname r=$r s=$s ... " | tee -a "$VERIFY_LOG"

        OUT=$(PIVOTER_COMPARE=1 OMP_NUM_THREADS=$THREADS timeout $TIMEOUT_SEC \
            "$BIN" "$graph" $r $s degen 2>&1)
        EXIT_CODE=$?

        if echo "$OUT" | grep -q "correctness verified\|Comparison passed"; then
            REF_TIME=$(extract_took "$OUT" "Reference")
            OPT_TIME=$(extract_took "$OUT" "Optimized")
            echo "PASS (ref=${REF_TIME}ms opt=${OPT_TIME}ms)" | tee -a "$VERIFY_LOG"
            VERIFY_PASS=$((VERIFY_PASS + 1))
        elif [ $EXIT_CODE -eq 124 ]; then
            echo "TIMEOUT" | tee -a "$VERIFY_LOG"
            VERIFY_SKIP=$((VERIFY_SKIP + 1))
        elif echo "$OUT" | grep -q "MISMATCH\|Comparison failed"; then
            echo "FAIL!" | tee -a "$VERIFY_LOG"
            echo "$OUT" >> "$VERIFY_LOG"
            VERIFY_FAIL=$((VERIFY_FAIL + 1))
        else
            echo "ERROR (exit=$EXIT_CODE)" | tee -a "$VERIFY_LOG"
            echo "$OUT" | tail -10 >> "$VERIFY_LOG"
            VERIFY_SKIP=$((VERIFY_SKIP + 1))
        fi
    done
done

echo "" | tee -a "$VERIFY_LOG"
echo "Verification: PASS=$VERIFY_PASS FAIL=$VERIFY_FAIL SKIP=$VERIFY_SKIP" | tee -a "$VERIFY_LOG"

if [ $VERIFY_FAIL -gt 0 ]; then
    echo "!!! CORRECTNESS FAILURES - ABORTING !!!" | tee -a "$VERIFY_LOG"
    exit 1
fi

# ============================================================
# Part 2: Performance benchmark (3 runs, median)
# ============================================================
PERF_LOG="$LOGDIR/performance_benchmark.log"
PERF_CSV="$LOGDIR/performance_benchmark.csv"

echo "" | tee "$PERF_LOG"
echo "============================================================" | tee -a "$PERF_LOG"
echo "  PART 2: Performance Benchmark" | tee -a "$PERF_LOG"
echo "============================================================" | tee -a "$PERF_LOG"

echo "graph,r,s,run1_ms,run2_ms,run3_ms,median_ms" > "$PERF_CSV"

for gi in "${!GRAPHS[@]}"; do
    graph="${GRAPHS[$gi]}"
    gname="${GRAPH_NAMES[$gi]}"
    [ ! -f "$graph" ] && continue

    for rs in "${RS_PAIRS[@]}"; do
        r=$(echo $rs | cut -d' ' -f1)
        s=$(echo $rs | cut -d' ' -f2)
        echo -n "$gname r=$r s=$s: " | tee -a "$PERF_LOG"

        TIMES=()
        for run in 1 2 3; do
            OUT=$(OMP_NUM_THREADS=$THREADS timeout $TIMEOUT_SEC \
                "$BIN" "$graph" $r $s degen 2>&1)
            if [ $? -eq 124 ]; then
                TIMES+=("TIMEOUT")
                echo -n "T " | tee -a "$PERF_LOG"
            else
                T=$(echo "$OUT" | grep "^time:" | awk '{print $2}')
                if [ -z "$T" ]; then
                    TIMES+=("ERR")
                    echo -n "E " | tee -a "$PERF_LOG"
                else
                    TIMES+=("$T")
                    echo -n "${T}ms " | tee -a "$PERF_LOG"
                fi
            fi
        done

        # Median
        NUMS=()
        for t in "${TIMES[@]}"; do
            [[ "$t" =~ ^[0-9]+$ ]] && NUMS+=("$t")
        done
        if [ ${#NUMS[@]} -ge 1 ]; then
            SORTED=($(printf '%s\n' "${NUMS[@]}" | sort -n))
            MEDIAN="${SORTED[$(( ${#SORTED[@]} / 2 ))]}"
        else
            MEDIAN="N/A"
        fi
        echo "-> median=${MEDIAN}ms" | tee -a "$PERF_LOG"
        echo "$gname,$r,$s,${TIMES[0]:-N/A},${TIMES[1]:-N/A},${TIMES[2]:-N/A},$MEDIAN" >> "$PERF_CSV"
    done
done

# ============================================================
# Part 3: Thread scaling
# ============================================================
SCALE_LOG="$LOGDIR/thread_scaling.log"
SCALE_CSV="$LOGDIR/thread_scaling.csv"

echo "" | tee "$SCALE_LOG"
echo "============================================================" | tee -a "$SCALE_LOG"
echo "  PART 3: Thread Scaling" | tee -a "$SCALE_LOG"
echo "============================================================" | tee -a "$SCALE_LOG"

THREAD_COUNTS="1 2 4 8 16 32 64"
SCALE_RS=("1 3" "2 3" "3 4" "4 5")

echo "graph,r,s,T1,T2,T4,T8,T16,T32,T64" > "$SCALE_CSV"

for gi in "${!GRAPHS[@]}"; do
    graph="${GRAPHS[$gi]}"
    gname="${GRAPH_NAMES[$gi]}"
    [ ! -f "$graph" ] && continue

    for rs in "${SCALE_RS[@]}"; do
        r=$(echo $rs | cut -d' ' -f1)
        s=$(echo $rs | cut -d' ' -f2)
        echo -n "$gname r=$r s=$s: " | tee -a "$SCALE_LOG"
        CSV_ROW="$gname,$r,$s"

        for T in $THREAD_COUNTS; do
            RESULT=$(OMP_NUM_THREADS=$T timeout $TIMEOUT_SEC \
                "$BIN" "$graph" $r $s degen 2>&1 | grep "^time:" | awk '{print $2}')
            if [ -z "$RESULT" ]; then
                RESULT="TIMEOUT"
                echo -n "T=$T:T " | tee -a "$SCALE_LOG"
            else
                echo -n "T=$T:${RESULT}ms " | tee -a "$SCALE_LOG"
            fi
            CSV_ROW="$CSV_ROW,$RESULT"
        done
        echo "" | tee -a "$SCALE_LOG"
        echo "$CSV_ROW" >> "$SCALE_CSV"
    done
done

# ============================================================
# Part 4: Optimized vs Reference timing
# ============================================================
REF_LOG="$LOGDIR/reference_comparison.log"
REF_CSV="$LOGDIR/reference_comparison.csv"

echo "" | tee "$REF_LOG"
echo "============================================================" | tee -a "$REF_LOG"
echo "  PART 4: Optimized vs Reference" | tee -a "$REF_LOG"
echo "============================================================" | tee -a "$REF_LOG"

echo "graph,r,s,ref_ms,opt_ms,speedup" > "$REF_CSV"

for gi in "${!GRAPHS[@]}"; do
    graph="${GRAPHS[$gi]}"
    gname="${GRAPH_NAMES[$gi]}"
    [ ! -f "$graph" ] && continue

    for rs in "${RS_PAIRS[@]}"; do
        r=$(echo $rs | cut -d' ' -f1)
        s=$(echo $rs | cut -d' ' -f2)
        echo -n "$gname r=$r s=$s: " | tee -a "$REF_LOG"

        OUT=$(PIVOTER_COMPARE=1 OMP_NUM_THREADS=$THREADS timeout $TIMEOUT_SEC \
            "$BIN" "$graph" $r $s degen 2>&1)

        if [ $? -eq 124 ]; then
            echo "TIMEOUT" | tee -a "$REF_LOG"
            echo "$gname,$r,$s,TIMEOUT,TIMEOUT,N/A" >> "$REF_CSV"
            continue
        fi

        REF_TIME=$(extract_took "$OUT" "Reference")
        OPT_TIME=$(extract_took "$OUT" "Optimized")

        if [ -n "$REF_TIME" ] && [ -n "$OPT_TIME" ]; then
            SPEEDUP=$(python3 -c "print(f'{${REF_TIME}/${OPT_TIME}:.2f}x')" 2>/dev/null || echo "N/A")
            echo "ref=${REF_TIME}ms opt=${OPT_TIME}ms speedup=${SPEEDUP}" | tee -a "$REF_LOG"
            echo "$gname,$r,$s,$REF_TIME,$OPT_TIME,$SPEEDUP" >> "$REF_CSV"
        else
            echo "parse error" | tee -a "$REF_LOG"
            echo "$OUT" | tail -5 >> "$REF_LOG"
            echo "$gname,$r,$s,ERR,ERR,N/A" >> "$REF_CSV"
        fi
    done
done

echo ""
echo "============================================================"
echo "  ALL EXPERIMENTS COMPLETE - $(date)"
echo "  Results: $LOGDIR"
echo "============================================================"
