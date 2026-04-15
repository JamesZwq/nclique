#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/.."

BIN="./build/bin/degeneracy_cliques"
TIMEOUT="${TIMEOUT:-3600}"
BASECSV="${BASECSV:-benchmark_all_results.csv}"
OUTCSV="${OUTCSV:-benchmark_recheck_results.csv}"
MERGEDCSV="${MERGEDCSV:-benchmark_all_results_merged.csv}"
LOGDIR="${LOGDIR:-bench_logs_recheck}"
DATADIR="${DATADIR:-/data/wenqianz}"
LOCKFILE="${LOCKFILE:-/tmp/bench_recheck_csv.lock}"
NPROC="${NPROC:-8}"

# Default mode:
#   134        rerun only ERROR(exit=134)
#   137        rerun only ERROR(exit=137)
#   all-error  rerun all ERROR(...)
#   timeout    rerun all TIMEOUT
MODE="${1:-134}"

run_one() {
  local graph="$1"
  local r="$2"
  local s="$3"
  local algo_label="$4"
  local env_var="$5"
  local logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_label}.log"

  if grep -q "^${graph},${r},${s},${algo_label}," "$OUTCSV" 2>/dev/null; then
    exit 0
  fi

  local result=""
  local exit_code=0
  result=$(env "${env_var}=1" timeout "${TIMEOUT}s" "$BIN" "graphs/${graph}.edges" "$r" "$s" 2>&1) || exit_code=$?
  echo "$result" > "$logfile"

  local took=""
  local status="OK"
  if [ "$exit_code" -eq 124 ]; then
    took="TIMEOUT"
    status="TIMEOUT"
  else
    took=$(echo "$result" | grep -oP 'NucleusCoreDecomposition took: \K[0-9.]+(?= ms)' || true)
    if [ -z "$took" ]; then
      took=$(echo "$result" | grep -oP 'ST_V\S+ took: \K[0-9.]+(?= ms)' || true)
    fi
    if [ -z "$took" ]; then
      took=$(echo "$result" | grep '^time:' | tail -1 | awk '{print $2}' || true)
    fi
    if [ -z "$took" ]; then
      took="ERROR"
      status="ERROR(exit=${exit_code})"
    fi
  fi

  local mem
  mem=$(echo "$result" | grep -oP 'Final Memory: \s*\K[0-9]+' || true)
  if [ -z "$mem" ]; then
    mem="N/A"
  fi

  (
    flock -x 200
    echo "${graph},${r},${s},${algo_label},${took},${mem},${status}" >> "$OUTCSV"
  ) 200>"$LOCKFILE"

  printf "    %-12s r=%s s=%-3s %10sms (%s) mem=%s kB\n" \
    "$algo_label" "$r" "$s" "$took" "$status" "$mem"
}

if [ "${1:-}" = "--run" ]; then
  run_one "$2" "$3" "$4" "$5" "$6"
  exit 0
fi

echo "============================================================="
echo "  Step 0: Incremental build"
echo "============================================================="
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -3
cmake --build build -j 12 --target degeneracy_cliques 2>&1 | tail -5
echo "Build done."
echo ""

echo "============================================================="
echo "  Step 1: Setup"
echo "============================================================="
mkdir -p graphs
for f in com-dblp.edges com-youtube.edges web-Stanford.edges web-Google.edges soc-pokec-relationships.edges web-it-2004.edges; do
  if [ ! -f "graphs/$f" ] && [ -f "$DATADIR/$f" ]; then
    ln -sf "$DATADIR/$f" "graphs/$f"
    echo "Linked graphs/$f -> $DATADIR/$f"
  fi
done

mkdir -p "$LOGDIR"
if [ ! -f "$OUTCSV" ] || ! head -1 "$OUTCSV" | grep -q 'graph,r,s'; then
  echo "graph,r,s,algorithm,time_ms,memory_kB,status" > "$OUTCSV"
fi

if [ ! -f "$BASECSV" ]; then
  echo "Base CSV not found: $BASECSV"
  exit 1
fi

JOBFILE=$(mktemp /tmp/bench_recheck_jobs.XXXXXX)

python3 - "$BASECSV" "$MODE" > "$JOBFILE" <<'PY'
import csv
import sys

basecsv, mode = sys.argv[1], sys.argv[2]

env_map = {
    "Ours_ST": "PIVOTER_RUN_ST",
    "REF_R1": "PIVOTER_NOOPT",
    "Ours_DCLP": "PIVOTER_RUN_R2_DCLP",
    "REF_R2": "PIVOTER_NOOPT",
    "Ours_V18": "PIVOTER_RUN_ST_V18",
    "Ours_V20": "PIVOTER_RUN_ST_V20",
    "REF_R3": "PIVOTER_RUN_REF",
}

want_timeout = mode == "timeout"
want_all_error = mode in {"all-error", "all_error", "all"}
want_codes = set()
if not want_timeout and not want_all_error:
    want_codes = {x.strip() for x in mode.split(",") if x.strip()}

seen = set()
with open(basecsv, newline="") as f:
    for row in csv.DictReader(f):
        status = row["status"]
        keep = False

        if want_timeout:
            keep = (status == "TIMEOUT")
        elif want_all_error:
            keep = status.startswith("ERROR(")
        elif status.startswith("ERROR(exit=") and status.endswith(")"):
            code = status[len("ERROR(exit="):-1]
            keep = code in want_codes

        if not keep:
            continue

        algo = row["algorithm"]
        env = env_map.get(algo)
        if env is None:
            continue

        key = (row["graph"], row["r"], row["s"], algo)
        if key in seen:
            continue
        seen.add(key)

        print("--run", row["graph"], row["r"], row["s"], algo, env)
PY

TOTAL=$(wc -l < "$JOBFILE" | tr -d ' ')
DONE=$(tail -n +2 "$OUTCSV" 2>/dev/null | wc -l | tr -d ' ')

echo "============================================================="
echo "  Step 2: Generate recheck jobs"
echo "============================================================="
echo "Mode: $MODE"
echo "Jobs to run: $TOTAL"
echo "Existing recheck rows: $DONE"
echo "Threads: $NPROC"
echo ""

if [ "$TOTAL" -eq 0 ]; then
  echo "No jobs matched."
else
  cat "$JOBFILE" | xargs -P "$NPROC" -L 1 bash "$(realpath "$0")"
fi

python3 - "$BASECSV" "$OUTCSV" "$MERGEDCSV" <<'PY'
import csv
import sys

basecsv, recheckcsv, mergedcsv = sys.argv[1], sys.argv[2], sys.argv[3]
fieldnames = ["graph", "r", "s", "algorithm", "time_ms", "memory_kB", "status"]

rows = {}
order = []

for path in [basecsv, recheckcsv]:
    try:
        f = open(path, newline="")
    except FileNotFoundError:
        continue
    with f:
        for row in csv.DictReader(f):
            key = (row["graph"], row["r"], row["s"], row["algorithm"])
            if key not in rows:
                order.append(key)
            rows[key] = row

with open(mergedcsv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for key in order:
        w.writerow(rows[key])

print(len(rows))
PY

rm -f "$JOBFILE" "$LOCKFILE"

echo ""
echo "============================================================="
echo "  DONE"
echo "============================================================="
echo "Recheck CSV: $OUTCSV"
echo "Merged CSV:  $MERGEDCSV"
echo "Logs:        $LOGDIR/"
echo "============================================================="
