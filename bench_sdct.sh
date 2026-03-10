#!/bin/bash
# 测试 SDCT 和 SDCT_Par 在不同线程数下的速度
#
# 用法: ./bench_sdct.sh <图文件目录> <最大线程数>
#
# 示例:
#   ./bench_sdct.sh graphs/ 8
#   ./bench_sdct.sh . 32

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <graph_dir> <max_threads>"
    echo "  graph_dir   : directory containing .edges files"
    echo "  max_threads : maximum thread count (doubles from 1 up to this)"
    exit 1
fi

GRAPH_DIR="$1"
MAX_THREADS="$2"
BIN="./build-ultra/bin/test_sdct_speed"
RESULT_FILE="bench_sdct_results_$(date +%Y%m%d_%H%M%S).txt"

# 检查二进制
if [ ! -f "$BIN" ]; then
    echo "[ERROR] Binary not found: $BIN"
    echo "  Please build first: cd build-ultra && make -j\$(nproc) test_sdct_speed"
    exit 1
fi

# 收集所有 .edges 文件
GRAPHS=()
while IFS= read -r f; do
    GRAPHS+=("$f")
done < <(find "$GRAPH_DIR" -maxdepth 1 -name '*.edges' | sort)

if [ ${#GRAPHS[@]} -eq 0 ]; then
    echo "[ERROR] No .edges files found in: $GRAPH_DIR"
    exit 1
fi

echo "=========================================="
echo "SDCT vs SDCT_Par Speed Benchmark"
echo "=========================================="
echo "Graph dir   : $GRAPH_DIR"
echo "Max threads : $MAX_THREADS"
echo "Graphs found: ${#GRAPHS[@]}"
echo "Binary      : $BIN"
echo "Results     : $RESULT_FILE"
echo ""

# 打印表头
HEADER=$(printf "%-40s  %-10s  %-10s  %-10s  %s" \
    "graph" "algorithm" "threads" "time_ms" "speedup")
echo "$HEADER" | tee "$RESULT_FILE"
echo "$(printf '%0.s-' {1..85})" | tee -a "$RESULT_FILE"

for GRAPH in "${GRAPHS[@]}"; do
    BASENAME=$(basename "$GRAPH")
    echo "" | tee -a "$RESULT_FILE"
    echo "--- $BASENAME ---" | tee -a "$RESULT_FILE"

    # 运行测试，解析每行输出
    "$BIN" "$GRAPH" "$MAX_THREADS" 2>/dev/null | while IFS= read -r line; do
        # 跳过 graph= 行，直接输出
        if [[ "$line" == graph=* ]]; then
            continue
        fi

        # 解析: SDCT threads=1 time_ms=123
        # 或:   SDCT_Par threads=4 time_ms=45 speedup=2.73
        if [[ "$line" =~ ^(SDCT[_A-Za-z]*)[[:space:]]threads=([0-9]+)[[:space:]]time_ms=([0-9]+)(.*) ]]; then
            ALG="${BASH_REMATCH[1]}"
            THR="${BASH_REMATCH[2]}"
            MS="${BASH_REMATCH[3]}"
            REST="${BASH_REMATCH[4]}"

            SPEEDUP="-"
            if [[ "$REST" =~ speedup=([0-9.]+) ]]; then
                SPEEDUP="${BASH_REMATCH[1]}x"
            fi

            ROW=$(printf "%-40s  %-10s  %-10s  %-10s  %s" \
                "$BASENAME" "$ALG" "$THR" "$MS" "$SPEEDUP")
            echo "$ROW" | tee -a "$RESULT_FILE"
        fi
    done
done

echo "" | tee -a "$RESULT_FILE"
echo "=========================================="  | tee -a "$RESULT_FILE"
echo "Done. Results saved to: $RESULT_FILE" | tee -a "$RESULT_FILE"
echo "=========================================="  | tee -a "$RESULT_FILE"
