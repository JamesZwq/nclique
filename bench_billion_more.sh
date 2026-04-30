#!/usr/bin/env bash
# Continue friendster billion-edge bench at s=5..30. Skips already-done rows.
# Stops early on memory pressure (>480GB peak) or empty CPI (no leaves).
# Run on tods2: cd ~/nclique && nohup bash bench_billion_more.sh > /tmp/billion_more.log 2>&1 &

set -u
ROOT="$HOME/nclique"
LINK="$ROOT/graphs/com-friendster.edges"
BIN="$ROOT/build/bin/degeneracy_cliques"
OUT_CSV="$ROOT/paper_data/bench_billion.csv"
LOG_DIR="$ROOT/bench_billion_logs"
mkdir -p "$LOG_DIR" "$ROOT/paper_data"
cd "$ROOT"

# Header if missing.
if [ ! -f "$OUT_CSV" ]; then
    echo "graph,s,algo,status,wall_sec,build_ms,peel_ms,sigma,mem_kB" > "$OUT_CSV"
fi

# Memory ceiling (kB). 480 GB.
MEM_CEILING_KB=503316480
# Upper bound on s. Practical: friendster has cliques but they thin out by s~25.
S_MAX=30

for s in $(seq 5 $S_MAX); do
    # Skip if already OK in CSV.
    if grep -q "^com-friendster,${s},V3,OK," "$OUT_CSV"; then
        echo "[$(date '+%F %T')] s=$s already OK, skipping"
        continue
    fi
    LOG="$LOG_DIR/friendster_s${s}_V3.log"
    echo "=== [$(date '+%F %T')] running s=$s ==="
    t0=$(date +%s)
    env PIVOTER_RUN_ST_V3=1 OMP_NUM_THREADS=24 \
        "$BIN" "$LINK" 1 "$s" 2>&1 | tee "$LOG" | tail -30
    rc=${PIPESTATUS[0]}
    t1=$(date +%s)
    wall=$((t1 - t0))

    # Extract stats from log directly (more reliable than regex on tee output).
    # Build phase: SDCT+callback or "Build took". Peel: ST_V3 r=1 (peel) took.
    build_ms=$(awk -F'[: ]+' '/SDCT\+callback took/ {print $4; exit}' "$LOG")
    peel_ms=$(awk -F'[: ]+' '/ST_V3 r=1 \(peel\) took/ {print $5; exit}' "$LOG")
    sigma=$(awk -F'=' '/COO entries=/ {gsub(/[^0-9]/, "", $NF); print $NF; exit}' "$LOG")
    mem=$(awk -F'[: \t]+' '/Final Memory:/ {print $3; exit}' "$LOG")

    status=OK
    if [ "$rc" -ne 0 ]; then status=ERR_$rc; fi

    echo "com-friendster,$s,V3,$status,$wall,${build_ms:-},${peel_ms:-},${sigma:-},${mem:-}" >> "$OUT_CSV"
    echo "[$(date '+%F %T')] s=$s done: $status wall=${wall}s sigma=${sigma:-?} mem=${mem:-?} kB"

    # Stop conditions.
    if [ "$status" != "OK" ]; then
        echo "[stop] s=$s failed; aborting"
        break
    fi
    # Empty CPI ⇒ no s-cliques exist, all higher s also empty.
    if [ -n "${sigma:-}" ] && [ "${sigma:-0}" -eq 0 ]; then
        echo "[stop] s=$s yielded sigma=0 (no $s-cliques); higher s vacuous"
        break
    fi
    # Memory pressure.
    if [ -n "${mem:-}" ] && [ "${mem:-0}" -gt "$MEM_CEILING_KB" ]; then
        echo "[stop] s=$s peak mem ${mem} kB exceeded ceiling ${MEM_CEILING_KB} kB"
        break
    fi
done

echo "=== ALL DONE ==="
cat "$OUT_CSV"
