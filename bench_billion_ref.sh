#!/usr/bin/env bash
# Run REF (the mutable-CPI baseline of NuclearCD) on com-friendster
# (1.806B edges) at s=2..S_MAX. Sister script to bench_billion.sh — same
# host, same machine, same /usr/bin/time -v wrapping, so the V3 vs REF
# rows can be diffed cell-for-cell.
#
# REF is selected by passing NO env flag — the binary's r=1 default is
# NCliqueVertexCoreDecomposition (the mutable BK tree implementation).
#
# Per the §7.3 finding, REF's peak RSS is up to 2.21× V3's, so high-s
# cells will OOM earlier than V3 did (V3 went 251 GB at s=4 on a 503 GB
# host with other users at ~160 GB).  Loop continues past OOM; sigma=0
# is the only "no more cliques exist" terminator.
#
# Usage:
#   ssh tods2; cd ~/nclique
#   nohup bash bench_billion_ref.sh > /tmp/billion_ref.log 2>&1 &
#
# Outputs: paper_data/friendster_billion/bench_billion_ref.csv
#          friendster_billion_ref_logs/

set -u
ROOT="$HOME/nclique"
GRAPH_DIR="/data/wenqianz"
GRAPH_GZ="$GRAPH_DIR/com-friendster.ungraph.txt.gz"
GRAPH="$GRAPH_DIR/com-friendster.ungraph.txt"
BIN="$ROOT/build/bin/degeneracy_cliques"
OUT_CSV="$ROOT/paper_data/friendster_billion/bench_billion_ref.csv"
LOG_DIR="$ROOT/friendster_billion_ref_logs"
TIME_BIN="/usr/bin/time"

S_MIN=${S_MIN:-2}
S_MAX=${S_MAX:-50}

mkdir -p "$LOG_DIR" "$ROOT/paper_data/friendster_billion"
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

# Stage 2: link into local graphs/.
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

# Stage 4: write CSV header if new — schema mirrors bench_billion.csv so
# the V3/REF rows can live in either file or be unioned trivially.
if [ ! -f "$OUT_CSV" ]; then
    echo "graph,s,algo,status,wall_sec,build_ms,peel_ms,sigma,mem_kB,time_max_rss_kB,time_user_sec,time_sys_sec,time_elapsed,time_pagefaults_major,time_pagefaults_minor" > "$OUT_CSV"
fi

# Stage 5: loop.
for s in $(seq "$S_MIN" "$S_MAX"); do
    if grep -qE "^com-friendster,${s},REF,(OK|ERR_137)," "$OUT_CSV"; then
        echo "[$(date '+%F %T')] s=$s already attempted (OK or OOM), skipping"
        continue
    fi

    LOG="$LOG_DIR/friendster_s${s}_REF.log"
    # NO env flag — that selects the REF default for r=1.  Still pin to 24
    # threads so any internal parallelism (counting, etc.) matches V3 runs.
    CMD_LINE="$TIME_BIN -v env OMP_NUM_THREADS=24 $BIN $LINK 1 $s"

    echo "=== [$(date '+%F %T')] running s=$s (REF) ==="
    echo "[CMD] $CMD_LINE" | tee "$LOG"

    t0=$(date +%s)
    "$TIME_BIN" -v env OMP_NUM_THREADS=24 "$BIN" "$LINK" 1 "$s" 2>&1 | tee -a "$LOG"
    rc=${PIPESTATUS[0]}
    t1=$(date +%s)
    wall=$((t1 - t0))

    # REF prints:
    #   "SDCT_Fused took: <ms> ms"            -> build phase
    #   "NucleusCoreDecomposition took: <ms> ms" -> peel phase (counting+peel)
    build_ms=$(awk -F'[: ]+' '/SDCT_Fused took/ {print $3; exit}' "$LOG")
    peel_ms=$(awk -F'[: ]+' '/NucleusCoreDecomposition took/ {print $3; exit}' "$LOG")
    sigma=$(awk -F'=' '/COO entries=/ {gsub(/[^0-9]/, "", $NF); print $NF; exit}' "$LOG")
    mem=$(awk -F'[: \t]+' '/Final Memory:/ {print $3; exit}' "$LOG")

    t_rss=$(awk -F': ' '/Maximum resident set size/ {print $2; exit}' "$LOG")
    t_user=$(awk -F': ' '/User time \(seconds\)/ {print $2; exit}' "$LOG")
    t_sys=$(awk -F': ' '/System time \(seconds\)/ {print $2; exit}' "$LOG")
    t_elapsed=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $NF; exit}' "$LOG")
    t_pfmaj=$(awk -F': ' '/Major \(requiring I\/O\) page faults/ {print $2; exit}' "$LOG")
    t_pfmin=$(awk -F': ' '/Minor \(reclaiming a frame\) page faults/ {print $2; exit}' "$LOG")

    status=OK
    if [ "$rc" -ne 0 ]; then status=ERR_$rc; fi

    echo "com-friendster,$s,REF,$status,$wall,${build_ms:-},${peel_ms:-},${sigma:-},${mem:-},${t_rss:-},${t_user:-},${t_sys:-},${t_elapsed:-},${t_pfmaj:-},${t_pfmin:-}" >> "$OUT_CSV"
    echo "[$(date '+%F %T')] s=$s done: $status wall=${wall}s build=${build_ms:-?}ms peel=${peel_ms:-?}ms maxRSS=${t_rss:-?}kB"

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
