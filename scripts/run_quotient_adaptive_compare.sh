#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <graph> <r> <s> [rounds] [tau_start] [tau_max] [lookahead] [auto_threshold_keys] [auto_raw_threshold]"
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT_DIR/build/bin/degeneracy_cliques"
GRAPH="$1"
R="$2"
S="$3"
ROUNDS="${4:-30}"
TAU_START="${5:-1}"
TAU_MAX="${6:-3}"
LOOKAHEAD="${7:-2}"
AUTO_THRESHOLD_KEYS="${8:-300000}"
AUTO_RAW_THRESHOLD="${9:-1000000}"

if [[ ! -x "$BIN" ]]; then
  echo "Binary not found: $BIN"
  echo "Build first with: cmake --build \"$ROOT_DIR/build\" -j8 --target degeneracy_cliques"
  exit 1
fi

if [[ ! -f "$GRAPH" ]]; then
  echo "Graph not found: $GRAPH"
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="$ROOT_DIR/benchmark_results/quotient_compare_${STAMP}"
mkdir -p "$OUT_DIR"

ADAPTIVE_LOG="$OUT_DIR/adaptive.log"
BUFFERED_LOG="$OUT_DIR/buffered.log"
AUTO_BUFFERED_LOG="$OUT_DIR/auto_buffered.log"
DYNAMIC_LOG="$OUT_DIR/dynamic.log"
SUMMARY_TXT="$OUT_DIR/summary.txt"

BANDED_AUTO_BUFFERED_LOG="$OUT_DIR/banded_auto_buffered.log"
BUCKETED_BANDED_AUTO_BUFFERED_LOG="$OUT_DIR/bucketed_banded_auto_buffered.log"

echo "[1/6] Running adaptive low-support prototype..."
PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
PIVOTER_QUOTIENT_LAB_ONLY=1 \
PIVOTER_QUOTIENT_ADAPTIVE_LOW=1 \
PIVOTER_QUOTIENT_ADAPTIVE_LOW_TAU_START="$TAU_START" \
PIVOTER_QUOTIENT_ADAPTIVE_LOW_TAU_MAX="$TAU_MAX" \
PIVOTER_QUOTIENT_MULTI_ROUND_MAX="$ROUNDS" \
"$BIN" "$GRAPH" "$R" "$S" > "$ADAPTIVE_LOG"

echo "[2/6] Running buffered adaptive low-support prototype..."
PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
PIVOTER_QUOTIENT_LAB_ONLY=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START="$TAU_START" \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX="$TAU_MAX" \
PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD="$LOOKAHEAD" \
PIVOTER_QUOTIENT_MULTI_ROUND_MAX="$ROUNDS" \
"$BIN" "$GRAPH" "$R" "$S" > "$BUFFERED_LOG"

echo "[3/6] Running auto-buffered adaptive low-support prototype..."
PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
PIVOTER_QUOTIENT_LAB_ONLY=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START="$TAU_START" \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX="$TAU_MAX" \
PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD="$LOOKAHEAD" \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_THRESHOLD_KEYS="$AUTO_THRESHOLD_KEYS" \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_RAW_THRESHOLD="$AUTO_RAW_THRESHOLD" \
PIVOTER_QUOTIENT_MULTI_ROUND_MAX="$ROUNDS" \
"$BIN" "$GRAPH" "$R" "$S" > "$AUTO_BUFFERED_LOG"

echo "[4/6] Running banded auto-buffered adaptive low-support prototype..."
PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
PIVOTER_QUOTIENT_LAB_ONLY=1 \
PIVOTER_QUOTIENT_BUFFERED_BANDED_LOW=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START="$TAU_START" \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX="$TAU_MAX" \
PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD="$LOOKAHEAD" \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_THRESHOLD_KEYS="$AUTO_THRESHOLD_KEYS" \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_RAW_THRESHOLD="$AUTO_RAW_THRESHOLD" \
PIVOTER_QUOTIENT_MULTI_ROUND_MAX="$ROUNDS" \
"$BIN" "$GRAPH" "$R" "$S" > "$BANDED_AUTO_BUFFERED_LOG"

echo "[5/6] Running bucketed banded auto-buffered adaptive low-support prototype..."
PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
PIVOTER_QUOTIENT_LAB_ONLY=1 \
PIVOTER_QUOTIENT_BUFFERED_BUCKETED_LOW=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO=1 \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_START="$TAU_START" \
PIVOTER_QUOTIENT_BUFFERED_LOW_TAU_MAX="$TAU_MAX" \
PIVOTER_QUOTIENT_BUFFERED_LOW_LOOKAHEAD="$LOOKAHEAD" \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_THRESHOLD_KEYS="$AUTO_THRESHOLD_KEYS" \
PIVOTER_QUOTIENT_BUFFERED_LOW_AUTO_RAW_THRESHOLD="$AUTO_RAW_THRESHOLD" \
PIVOTER_QUOTIENT_MULTI_ROUND_MAX="$ROUNDS" \
"$BIN" "$GRAPH" "$R" "$S" > "$BUCKETED_BANDED_AUTO_BUFFERED_LOG"

echo "[6/6] Running full dynamic sparse prototype..."
PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
PIVOTER_QUOTIENT_LAB_ONLY=1 \
PIVOTER_QUOTIENT_MULTI_ROUND_DYNAMIC=1 \
PIVOTER_QUOTIENT_MULTI_ROUND_MAX="$ROUNDS" \
"$BIN" "$GRAPH" "$R" "$S" > "$DYNAMIC_LOG"

{
  echo "Graph: $GRAPH"
  echo "r=$R s=$S rounds=$ROUNDS tau_start=$TAU_START tau_max=$TAU_MAX lookahead=$LOOKAHEAD auto_threshold_keys=$AUTO_THRESHOLD_KEYS auto_raw_threshold=$AUTO_RAW_THRESHOLD"
  echo
  echo "=== Adaptive summary ==="
  sed -n '/^========== Quotient Adaptive Low-Support ==========/,$p' "$ADAPTIVE_LOG" \
    | rg "Tau start/max|Tau phases used|Completed rounds|Total exact leaves|Total candidate keys|Total spawned subleaves|Max observed round min|Total rebuild time|Final tau|Prototype time|Phase tau=" || true
  echo
  echo "=== Buffered summary ==="
  sed -n '/^========== Quotient Buffered Adaptive Low-Support ==========/,$p' "$BUFFERED_LOG" \
    | rg "Tau start/max|Lookahead|Auto lookahead|Rebuild count|Completed rounds|Total exact leaves|Total candidate keys|Total spawned subleaves|Max observed round min|Total rebuild time|Final tau|Prototype time|Phase tau=" || true
  echo
  echo "=== Auto-buffered summary ==="
  sed -n '/^========== Quotient Buffered Adaptive Low-Support ==========/,$p' "$AUTO_BUFFERED_LOG" \
    | rg "Tau start/max|Lookahead|Auto lookahead|Auto threshold keys|Auto raw threshold|Rebuild count|Completed rounds|Total exact leaves|Total candidate keys|Total spawned subleaves|Max observed round min|Total rebuild time|Final tau|Prototype time|Phase tau=" || true
  echo
  echo "=== Banded auto-buffered summary ==="
  sed -n '/^======= Quotient Banded Buffered Adaptive Low-Support =======/,$p' "$BANDED_AUTO_BUFFERED_LOG" \
    | rg "Tau start/max|Lookahead|Auto lookahead|Auto threshold keys|Auto raw threshold|Rebuild count|Completed rounds|Total exact leaves|Total candidate keys|Total spawned subleaves|Max observed round min|Total rebuild time|Final tau|Prototype time|Phase tau=" || true
  echo
  echo "=== Bucketed banded auto-buffered summary ==="
  sed -n '/^==== Quotient Bucketed Banded Buffered Adaptive Low-Support ====/,$p' "$BUCKETED_BANDED_AUTO_BUFFERED_LOG" \
    | rg "Tau start/max|Lookahead|Auto lookahead|Auto threshold keys|Auto raw threshold|Rebuild count|Completed rounds|Total exact leaves|Total candidate keys|Total spawned subleaves|Max observed round min|Total rebuild time|Final tau|Prototype time|Phase tau=" || true
  echo
  echo "=== Dynamic summary ==="
  sed -n '/^========= Quotient Multi-Round Dynamic Prototype =========$/,$p' "$DYNAMIC_LOG" \
    | rg "Completed rounds|First core>1 round|First split round|Rounds with split|Dead leaves so far|Split leaves so far|Single-child leaves|Multi-child leaves|Spawned leaves so far|Max spawned in one round|Max observed core|Remaining sparse states|Prototype time" || true
} | tee "$SUMMARY_TXT"

echo
echo "Logs written to:"
echo "  $ADAPTIVE_LOG"
echo "  $BUFFERED_LOG"
echo "  $AUTO_BUFFERED_LOG"
echo "  $BANDED_AUTO_BUFFERED_LOG"
echo "  $BUCKETED_BANDED_AUTO_BUFFERED_LOG"
echo "  $DYNAMIC_LOG"
echo "  $SUMMARY_TXT"
