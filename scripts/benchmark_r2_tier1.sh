#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

BIN="${BIN:-./build/bin/degeneracy_cliques}"
R="${R:-2}"
S="${S:-3}"
REPEATS="${REPEATS:-3}"
GRAPHS="${GRAPHS:-com-youtube web-Stanford web-it-2004}"
OUTDIR="${OUTDIR:-benchmark_results/r2_tier1_$(date +%Y%m%d_%H%M%S)}"
CSV="$OUTDIR/results.csv"

mkdir -p "$OUTDIR"
echo "graph,algo,run,time_ms,memory_kB,status" > "$CSV"

run_one() {
  local graph="$1"
  local algo="$2"
  local run_id="$3"
  local logfile="$OUTDIR/${graph}_${algo}_run${run_id}.log"
  local result=""
  local status="OK"
  local exit_code=0

  if [[ "$algo" == "DCLP" ]]; then
    result=$(env PIVOTER_RUN_R2_DCLP=1 "$BIN" "graphs/${graph}.edges" "$R" "$S" 2>&1) || exit_code=$?
  elif [[ "$algo" == "HYBRID_LAB" ]]; then
    result=$(env PIVOTER_RUN_R2_HYBRID_LAB=1 "$BIN" "graphs/${graph}.edges" "$R" "$S" 2>&1) || exit_code=$?
  else
    echo "unknown algo: $algo" >&2
    exit 1
  fi

  echo "$result" > "$logfile"

  local parsed took mem
  parsed=$(python3 - "$logfile" <<'PY'
import pathlib, re, sys

text = pathlib.Path(sys.argv[1]).read_text()
took = re.findall(r'NucleusCoreDecomposition took:\s*([0-9.]+)\s*ms', text)
mem = re.findall(r'Final Memory:\s*([0-9]+)', text)
print((took[-1] if took else "") + "\t" + (mem[-1] if mem else ""))
PY
)
  took="${parsed%%$'\t'*}"
  mem="${parsed#*$'\t'}"
  if [[ "$mem" == "$parsed" ]]; then
    mem=""
  fi

  if [[ $exit_code -ne 0 ]]; then
    status="ERROR(exit=$exit_code)"
  fi
  if [[ -z "${took}" ]]; then
    took="ERROR"
    status="PARSE_FAIL"
  fi
  if [[ -z "${mem}" ]]; then
    mem="N/A"
  fi

  echo "$graph,$algo,$run_id,$took,$mem,$status" >> "$CSV"
  printf "%-14s %-11s run=%s time=%sms mem=%s status=%s\n" "$graph" "$algo" "$run_id" "$took" "$mem" "$status"
}

echo "Tier-1 benchmark"
echo "  graphs: $GRAPHS"
echo "  repeats: $REPEATS"
echo "  output: $OUTDIR"
echo ""

for graph in $GRAPHS; do
  for run_id in $(seq 1 "$REPEATS"); do
    if (( run_id % 2 == 1 )); then
      run_one "$graph" "DCLP" "$run_id"
      run_one "$graph" "HYBRID_LAB" "$run_id"
    else
      run_one "$graph" "HYBRID_LAB" "$run_id"
      run_one "$graph" "DCLP" "$run_id"
    fi
  done
done

echo ""
echo "Summary"
python3 - "$CSV" <<'PY'
import csv, math, statistics, sys
from collections import defaultdict

path = sys.argv[1]
rows = []
with open(path) as f:
    for row in csv.DictReader(f):
        if row["status"] == "OK":
            rows.append(row)

bucket = defaultdict(list)
for row in rows:
    bucket[(row["graph"], row["algo"])].append(float(row["time_ms"]))

medians = {}
for (graph, algo), vals in sorted(bucket.items()):
    vals = sorted(vals)
    median = statistics.median(vals)
    mean = statistics.fmean(vals)
    medians[(graph, algo)] = median
    print(f"{graph:14s} {algo:11s} n={len(vals)} min={vals[0]:.3f} median={median:.3f} mean={mean:.3f} max={vals[-1]:.3f}")

ratios = []
for graph in sorted({graph for graph, _ in medians}):
    dclp = medians.get((graph, "DCLP"))
    hybrid = medians.get((graph, "HYBRID_LAB"))
    if dclp is None or hybrid is None:
        continue
    ratio = hybrid / dclp
    ratios.append(ratio)
    print(f"ratio {graph:14s} HYBRID_LAB/DCLP = {ratio:.4f}x")

if ratios:
    geom = math.exp(sum(math.log(r) for r in ratios) / len(ratios))
    print(f"geom_mean_ratio HYBRID_LAB/DCLP = {geom:.4f}x")
PY

echo ""
echo "CSV: $CSV"
