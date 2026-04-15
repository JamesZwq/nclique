#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT_DIR/build/bin/degeneracy_cliques"
OUT_CSV="${1:-$ROOT_DIR/quotient_state_survey.csv}"

if [[ ! -x "$BIN" ]]; then
  echo "Binary not found: $BIN" >&2
  exit 1
fi

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

cat > "$OUT_CSV" <<'EOF'
graph,r,s,sdct_ms,leaf_count,total_clean,total_refined,total_delta,total_touch_pct,median_clean,p90_clean,median_refined,p90_refined,median_delta,p90_delta,median_clean_mult,p90_clean_mult,median_delta_mult,p90_delta_mult,median_touch_pct,p90_touch_pct,proto_mem_mb,avg_clean_per_leaf,avg_delta_per_leaf
EOF

declare -a CASES=(
  "graphs/com-youtube.edges 3 4"
  "graphs/com-youtube.edges 4 5"
  "graphs/web-Stanford.edges 4 5"
  "graphs/web-it-2004.edges 3 4"
  "graphs/web-it-2004.edges 4 5"
)

for item in "${CASES[@]}"; do
  read -r graph r s <<<"$item"
  base="$(basename "$graph" .edges)"
  log="$TMP_DIR/${base}_r${r}_s${s}.log"

  echo "[survey] $graph r=$r s=$s"
  env \
    PIVOTER_RUN_ST_QUOTIENT_LAB=1 \
    PIVOTER_QUOTIENT_LAB_ONLY=1 \
    PIVOTER_QUOTIENT_BUILD_STATE=1 \
    "$BIN" "$ROOT_DIR/$graph" "$r" "$s" > "$log" 2>&1

  python3 - "$log" "$graph" "$r" "$s" >> "$OUT_CSV" <<'PY'
import re
import sys

log_path, graph, r, s = sys.argv[1:]
text = open(log_path, "r", encoding="utf-8", errors="ignore").read()

def grab(pattern, default=""):
    m = re.search(pattern, text)
    return m.group(1) if m else default

def grab_pair(pattern):
    m = re.search(pattern, text)
    if not m:
        return ("", "")
    return (m.group(1), m.group(2))

sdct = grab(r"SDCT_Fused took:\s*([0-9.e+-]+)\s*ms")
leaf_count = grab(r"SDCT leaf count \(stored\):\s*([0-9]+)")
total_clean = grab(r"Total clean quotient compression:\s*([0-9.e+-]+)x")
total_refined = grab(r"Total one-step refined compression:\s*([0-9.e+-]+)x")
total_delta = grab(r"Total one-removed delta compression:\s*([0-9.e+-]+)x")
total_touch = grab(r"Total one-removed touch fraction:\s*([0-9.e+-]+)%")
median_clean, p90_clean = grab_pair(r"Per-leaf compression median/p90/max:\s*([0-9.e+-]+)x\s*/\s*([0-9.e+-]+)x")
median_refined, p90_refined = grab_pair(r"Per-leaf one-step refined median/p90/max:\s*([0-9.e+-]+)x\s*/\s*([0-9.e+-]+)x")
median_delta, p90_delta = grab_pair(r"Per-leaf one-removed delta median/p90/max:\s*([0-9.e+-]+)x\s*/\s*([0-9.e+-]+)x")
median_clean_mult, p90_clean_mult = grab_pair(r"Largest clean class multiplicity median/p90/max:\s*([0-9.e+-]+)\s*/\s*([0-9.e+-]+)\s*/")
median_delta_mult, p90_delta_mult = grab_pair(r"Largest delta class multiplicity median/p90/max:\s*([0-9.e+-]+)\s*/\s*([0-9.e+-]+)\s*/")
median_touch, p90_touch = grab_pair(r"One-removed touch fraction median/p90/max:\s*([0-9.e+-]+)%\s*/\s*([0-9.e+-]+)%\s*/")
proto_mem = grab(r"Prototype state memory:\s*([0-9]+)\s*MB")
avg_clean = grab(r"Avg clean entries/leaf:\s*([0-9.e+-]+)")
avg_delta = grab(r"Avg delta entries/leaf:\s*([0-9.e+-]+)")

print(",".join([
    graph, r, s, sdct, leaf_count, total_clean, total_refined, total_delta, total_touch,
    median_clean, p90_clean, median_refined, p90_refined, median_delta, p90_delta,
    median_clean_mult, p90_clean_mult, median_delta_mult, p90_delta_mult,
    median_touch, p90_touch,
    proto_mem, avg_clean, avg_delta
]))
PY
done

echo "Wrote $OUT_CSV"
