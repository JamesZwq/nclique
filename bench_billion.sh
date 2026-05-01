#!/usr/bin/env bash
# Run V3 on friendster (1.8B edges) at s=2..S_MAX. No timeout.
# Captures /usr/bin/time -v stats (wall, user, sys, peak RSS, page faults)
# and dumps per-vertex core values via PIVOTER_DUMP_CORE for archival.
# Skips s values already recorded as OK in the CSV (resume support).
#
# Usage:  ssh tods2; cd ~/nclique; nohup bash bench_billion.sh > /tmp/billion.log 2>&1 &
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
TIME_BIN="/usr/bin/time"

# Practical upper bound; bench stops early on ERR / Σ=0 / mem ceiling.
S_MIN=${S_MIN:-2}
S_MAX=${S_MAX:-50}
# 480 GB ceiling (host has 503 GB).
MEM_CEILING_KB=${MEM_CEILING_KB:-503316480}

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

# Stage 2: link into local graphs/ dir.
mkdir -p "$ROOT/graphs"
LINK="$ROOT/graphs/com-friendster.edges"
[ -e "$LINK" ] || ln -s "$GRAPH" "$LINK"

# Stage 3: rebuild.
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build "$ROOT/build" -j 12 --target degeneracy_cliques 2>&1 | tail -3

if [ ! -x "$TIME_BIN" ]; then
    echo "ERROR: $TIME_BIN not found (need GNU time, not bash builtin)"
    exit 1
fi

# Stage 4: write CSV header if new. Extra columns for /usr/bin/time -v.
if [ ! -f "$OUT_CSV" ]; then
    echo "graph,s,algo,status,wall_sec,build_ms,peel_ms,sigma,mem_kB,time_max_rss_kB,time_user_sec,time_sys_sec,time_elapsed,time_pagefaults_major,time_pagefaults_minor" > "$OUT_CSV"
fi

# Stage 5: loop.
for s in $(seq "$S_MIN" "$S_MAX"); do
    # Skip if already attempted (OK or OOM-killed by signal 137).
    # Other ERRs may be transient (e.g., interrupted) and worth retrying.
    if grep -qE "^com-friendster,${s},V3,(OK|ERR_137)," "$OUT_CSV"; then
        echo "[$(date '+%F %T')] s=$s already attempted (OK or OOM), skipping"
        continue
    fi

    LOG="$LOG_DIR/friendster_s${s}_V3.log"
    DUMP="$LOG_DIR/friendster_s${s}_V3.cores"
    CMD_LINE="$TIME_BIN -v env PIVOTER_RUN_ST_V3=1 OMP_NUM_THREADS=24 PIVOTER_DUMP_CORE=$DUMP $BIN $LINK 1 $s"

    echo "=== [$(date '+%F %T')] running s=$s ==="
    echo "[CMD] $CMD_LINE" | tee "$LOG"

    t0=$(date +%s)
    "$TIME_BIN" -v env PIVOTER_RUN_ST_V3=1 OMP_NUM_THREADS=24 PIVOTER_DUMP_CORE="$DUMP" \
        "$BIN" "$LINK" 1 "$s" 2>&1 | tee -a "$LOG"
    rc=${PIPESTATUS[0]}
    t1=$(date +%s)
    wall=$((t1 - t0))

    # Extract from log.
    build_ms=$(awk -F'[: ]+' '/SDCT\+callback took/ {print $4; exit}' "$LOG")
    peel_ms=$(awk -F'[: ]+' '/ST_V3 r=1 \(peel\) took/ {print $5; exit}' "$LOG")
    sigma=$(awk -F'=' '/COO entries=/ {gsub(/[^0-9]/, "", $NF); print $NF; exit}' "$LOG")
    mem=$(awk -F'[: \t]+' '/Final Memory:/ {print $3; exit}' "$LOG")

    # /usr/bin/time -v fields. Match against the verbose output's labels.
    t_rss=$(awk -F': ' '/Maximum resident set size/ {print $2; exit}' "$LOG")
    t_user=$(awk -F': ' '/User time \(seconds\)/ {print $2; exit}' "$LOG")
    t_sys=$(awk -F': ' '/System time \(seconds\)/ {print $2; exit}' "$LOG")
    t_elapsed=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $NF; exit}' "$LOG")
    t_pfmaj=$(awk -F': ' '/Major \(requiring I\/O\) page faults/ {print $2; exit}' "$LOG")
    t_pfmin=$(awk -F': ' '/Minor \(reclaiming a frame\) page faults/ {print $2; exit}' "$LOG")

    status=OK
    if [ "$rc" -ne 0 ]; then status=ERR_$rc; fi

    echo "com-friendster,$s,V3,$status,$wall,${build_ms:-},${peel_ms:-},${sigma:-},${mem:-},${t_rss:-},${t_user:-},${t_sys:-},${t_elapsed:-},${t_pfmaj:-},${t_pfmin:-}" >> "$OUT_CSV"
    echo "[$(date '+%F %T')] s=$s done: $status wall=${wall}s build=${build_ms:-?}ms peel=${peel_ms:-?}ms maxRSS=${t_rss:-?}kB"

    # Stop only if we have proven there are no more $s$-cliques (sigma=0):
    # at that point all higher s are vacuous. ERR (incl. OOM) does NOT stop the
    # loop --- as s grows past the max-clique size, Sigma can shrink again, so
    # higher s may fit even if a middle s did not.
    if [ -n "${sigma:-}" ] && [ "${sigma}" -eq 0 ]; then
        echo "[stop] s=$s yielded sigma=0; higher s vacuous"
        break
    fi
    if [ "$status" != "OK" ]; then
        echo "[continue] s=$s failed (rc=$rc); trying next s (CPI may shrink past the peak)"
    fi
done

echo "=== ALL DONE ==="
cat "$OUT_CSV"
