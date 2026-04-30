#!/usr/bin/env bash
# Run V3 on friendster (1.8B edges) to demonstrate billion-edge scalability.
# Targets s=2 (k-core, baseline), s=3 (triangles), and s=4 if memory permits.
#
# Usage:  ssh tods2; cd ~/nclique; bash bench_billion.sh
#
# Outputs to paper_data/bench_billion.csv and bench_billion_logs/.

set -u
ROOT="$HOME/nclique"
GRAPH_DIR="/data/wenqianz"
GRAPH_GZ="$GRAPH_DIR/com-friendster.ungraph.txt.gz"
GRAPH="$GRAPH_DIR/com-friendster.ungraph.txt"
BIN="$ROOT/build/bin/degeneracy_cliques"
OUT_CSV="$ROOT/paper_data/bench_billion.csv"
LOG_DIR="$ROOT/bench_billion_logs"

mkdir -p "$LOG_DIR" "$ROOT/paper_data"
cd "$ROOT"

# Stage 1: decompress if needed.
if [ ! -f "$GRAPH" ] && [ -f "$GRAPH_GZ" ]; then
    echo "[$(date '+%F %T')] decompressing $GRAPH_GZ ..."
    gunzip -k "$GRAPH_GZ"
fi
if [ ! -f "$GRAPH" ]; then
    echo "ERROR: $GRAPH not found and no .gz to decompress"
    exit 1
fi

ls -lh "$GRAPH"

# Stage 2: link into local graphs/ dir so V3 can find it.
mkdir -p "$ROOT/graphs"
LINK="$ROOT/graphs/com-friendster.edges"
[ -e "$LINK" ] || ln -s "$GRAPH" "$LINK"

# Stage 3: rebuild binary so it's fresh.
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build "$ROOT/build" -j 12 --target degeneracy_cliques 2>&1 | tail -3

# Stage 4: write CSV header if new.
if [ ! -f "$OUT_CSV" ]; then
    echo "graph,s,algo,status,wall_sec,build_ms,peel_ms,sigma,mem_kB" > "$OUT_CSV"
fi

# Stage 5: run V3 at s=2,3,4. No timeout — run to completion.
for s in 2 3 4; do
    LOG="$LOG_DIR/friendster_s${s}_V3.log"
    echo "=== [$(date '+%F %T')] running s=$s ==="
    t0=$(date +%s)
    env PIVOTER_RUN_ST_V3=1 OMP_NUM_THREADS=24 \
        "$BIN" "$LINK" 1 "$s" 2>&1 | tee "$LOG" | tail -30
    rc=$?
    t1=$(date +%s)
    wall=$((t1 - t0))

    # Extract stats
    build_ms=$(grep -oE "ST_V3 Build took: [0-9.]+ ms" "$LOG" | head -1 | grep -oE "[0-9.]+" | head -1)
    peel_ms=$(grep -oE "ST_V3 r=1 \(peel\) took: [0-9.]+ ms" "$LOG" | head -1 | grep -oE "[0-9.]+" | head -1)
    sigma=$(grep -oE "COO entries=[0-9]+" "$LOG" | head -1 | grep -oE "[0-9]+")
    mem=$(grep -oE "Final Memory:[[:space:]]+[0-9.]+ kB" "$LOG" | head -1 | grep -oE "[0-9.]+")
    status=OK
    if [ "$rc" -ne 0 ]; then status=ERR_$rc; fi
    echo "com-friendster,$s,V3,$status,$wall,${build_ms:-},${peel_ms:-},${sigma:-},${mem:-}" >> "$OUT_CSV"
    echo "[$(date '+%F %T')] s=$s done: $status wall=${wall}s build=${build_ms}ms peel=${peel_ms}ms"
    if [ "$status" != "OK" ]; then break; fi
done

echo "=== ALL DONE ==="
cat "$OUT_CSV"
