#!/bin/bash
# Benchmark script: test r=2 s=3 and r=1 s=3 at various thread counts
# Usage: ./benchmark_threads.sh <graph_dir> [binary_path]
#   e.g.: ./benchmark_threads.sh /path/to/graphs ./build/bin/degeneracy_cliques

GRAPH_DIR="${1:-.}"
BIN="${2:-./build/bin/degeneracy_cliques}"

if [ ! -f "$BIN" ]; then
    echo "Error: binary not found at $BIN"
    echo "Build first: mkdir -p build && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"
    exit 1
fi

THREADS="1 2 4 8 16 32 64"
GRAPHS="email-Eu-core email-Enron com-dblp web-Stanford com-youtube soc-pokec-relationships"

echo "============================================================"
echo "  Thread Scaling Benchmark"
echo "  Binary: $BIN"
echo "  Graph dir: $GRAPH_DIR"
echo "  Date: $(date)"
echo "  CPU: $(uname -p) $(nproc 2>/dev/null || sysctl -n hw.ncpu) cores"
echo "============================================================"
echo ""

# r=2 s=3
echo "==================== r=2 s=3 ===================="
printf "%-30s" "Graph"
for T in $THREADS; do
    printf "%8s" "T=$T"
done
echo ""
printf "%-30s" "-----"
for T in $THREADS; do
    printf "%8s" "-----"
done
echo ""

for G in $GRAPHS; do
    GFILE="$GRAPH_DIR/$G.edges"
    if [ ! -f "$GFILE" ]; then
        continue
    fi
    printf "%-30s" "$G"
    for T in $THREADS; do
        RESULT=$(OMP_NUM_THREADS=$T "$BIN" "$GFILE" 2 3 degen 2>&1 | grep "^time:" | awk '{print $2}')
        if [ -z "$RESULT" ]; then
            printf "%8s" "ERR"
        else
            printf "%7s%s" "$RESULT" "ms"
        fi
    done
    echo ""
done

echo ""

# r=1 s=3
echo "==================== r=1 s=3 ===================="
printf "%-30s" "Graph"
for T in $THREADS; do
    printf "%8s" "T=$T"
done
echo ""
printf "%-30s" "-----"
for T in $THREADS; do
    printf "%8s" "-----"
done
echo ""

for G in $GRAPHS; do
    GFILE="$GRAPH_DIR/$G.edges"
    if [ ! -f "$GFILE" ]; then
        continue
    fi
    printf "%-30s" "$G"
    for T in $THREADS; do
        RESULT=$(OMP_NUM_THREADS=$T "$BIN" "$GFILE" 1 3 degen 2>&1 | grep "^time:" | awk '{print $2}')
        if [ -z "$RESULT" ]; then
            printf "%8s" "ERR"
        else
            printf "%7s%s" "$RESULT" "ms"
        fi
    done
    echo ""
done

echo ""

# r=2 s=4
echo "==================== r=2 s=4 ===================="
printf "%-30s" "Graph"
for T in $THREADS; do
    printf "%8s" "T=$T"
done
echo ""
printf "%-30s" "-----"
for T in $THREADS; do
    printf "%8s" "-----"
done
echo ""

for G in $GRAPHS; do
    GFILE="$GRAPH_DIR/$G.edges"
    if [ ! -f "$GFILE" ]; then
        continue
    fi
    printf "%-30s" "$G"
    for T in $THREADS; do
        RESULT=$(OMP_NUM_THREADS=$T "$BIN" "$GFILE" 2 4 degen 2>&1 | grep "^time:" | awk '{print $2}')
        if [ -z "$RESULT" ]; then
            printf "%8s" "ERR"
        else
            printf "%7s%s" "$RESULT" "ms"
        fi
    done
    echo ""
done

echo ""
echo "Done."
