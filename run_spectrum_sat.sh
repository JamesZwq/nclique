#!/usr/bin/env bash
# run_spectrum_sat.sh — drive bench_spectrum_sat over a graph list, capturing
# the full stdout+stderr per graph (peak RSS via /usr/bin/time) and collecting
# the SAT_SUMMARY lines into a CSV.
#
# Usage: ./run_spectrum_sat.sh <out_dir> <graph1.edges> [graph2.edges ...]
#   env: SMAX=<N> to cap s (default: binary's 400 cap); SORT=degen|degenR|default
#
# Build first:
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j 12 \
#         --target bench_spectrum_sat
set -u
BIN=./build/bin/bench_spectrum_sat
OUT="${1:?usage: run_spectrum_sat.sh <out_dir> <graphs...>}"; shift
mkdir -p "$OUT"
CSV="$OUT/spectrum_sat.csv"
SMAX_ARG=(); [ -n "${SMAX:-}" ] && SMAX_ARG=(--smax "$SMAX")
SORT_ARG=(--sort "${SORT:-degen}")

# macOS uses `time -l`; Linux servers use `time -v`. Detect.
if /usr/bin/time -l true 2>/dev/null; then TIMEFLAG=-l; else TIMEFLAG=-v; fi

echo "graph,n,m,maxdeg,omega,omega_trunc,n_active,saturated,sat_frac,nonzero_deltas,traj_len_sum,naive_store,anchordelta_store,comp_ratio,n_active_exact,sat_frac_exact,comp_ratio_exact,law_viol,nest_viol,fuzzy,total_ms" > "$CSV"

for f in "$@"; do
  [ -f "$f" ] || { echo "MISSING $f" >&2; continue; }
  base=$(basename "$f" .edges)
  log="$OUT/$base.log"
  echo ">>> $base"
  /usr/bin/time $TIMEFLAG "$BIN" "$f" "${SMAX_ARG[@]}" "${SORT_ARG[@]}" >"$log" 2>&1
  line=$(grep -m1 '^SAT_SUMMARY' "$log")
  [ -z "$line" ] && { echo "  (no SAT_SUMMARY — see $log)"; continue; }
  # parse key=val tokens into a CSV row, preserving column order
  echo "$line" | awk '{
    for (i=2;i<=NF;i++){ split($i,a,"="); v[a[1]]=a[2] }
    print v["graph"]","v["n"]","v["m"]","v["maxdeg"]","v["omega"]","v["omega_trunc"]","\
          v["n_active"]","v["saturated"]","v["sat_frac"]","v["nonzero_deltas"]","\
          v["traj_len_sum"]","v["naive_store"]","v["anchordelta_store"]","v["comp_ratio"]","\
          v["n_active_exact"]","v["sat_frac_exact"]","v["comp_ratio_exact"]","\
          v["law_viol"]","v["nest_viol"]","v["fuzzy"]","v["total_ms"]
  }' >> "$CSV"
  echo "  $line"
done
echo "CSV -> $CSV"
