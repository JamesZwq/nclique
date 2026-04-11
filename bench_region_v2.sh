#!/bin/bash
# =============================================================
# Region Tuple V2 Benchmark
# Compare: REF (Correct), ST (CPI Peeling), Region V2 (Tuple Peeling)
#
# Metrics: total time, phase breakdown, peak memory
# Output: CSV + human-readable table
#
# Usage:
#   bash bench_region_v2.sh               # run all
#   bash bench_region_v2.sh --quick       # quick subset only
#   bash bench_region_v2.sh --verify      # also check correctness
# =============================================================

set -uo pipefail
cd "$(dirname "$0")"

BIN="./build/bin/degeneracy_cliques"
TIMEOUT=3600
OUTDIR="benchmark_results/region_v2_$(date +%Y%m%d_%H%M%S)"
OUTCSV="$OUTDIR/results.csv"
LOGDIR="$OUTDIR/logs"
QUICK=0
VERIFY=0

for arg in "$@"; do
  case "$arg" in
    --quick) QUICK=1 ;;
    --verify) VERIFY=1 ;;
  esac
done

# ============ Build ============
echo "============================================================="
echo "  Region Tuple V2 Benchmark"
echo "============================================================="
echo ""
echo "Building..."
mkdir -p build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release 2>&1 | tail -1
cmake --build build -j 12 --target degeneracy_cliques 2>&1 | tail -1
echo "Build done."
echo ""

mkdir -p "$LOGDIR"
echo "graph,r,s,algorithm,total_ms,enum_ms,build_ms,peel_ms,peak_mem_kB,status,max_core,r_cliques,r_tuples,s_tuples" > "$OUTCSV"

# ============ Datasets ============
# Format: graph:omega (max clique size)
if [ $QUICK -eq 1 ]; then
  DATASETS=(
    "dblp-core30:114"
    "email-Eu-core:25"
    "com-dblp:113"
  )
  RS_COMBOS="3,4 3,5 4,5"
else
  DATASETS=(
    "dblp-core30:114"
    "email-Eu-core:25"
    "email-Enron:20"
    "com-dblp:113"
    "com-youtube:17"
    "web-Stanford:61"
    "web-it-2004:432"
    "soc-pokec-relationships:47"
  )
  RS_COMBOS="3,4 3,5 3,6 4,5 4,6 5,6"
fi

# ============ Algorithms ============
# name:env_var
ALGOS=(
  "REF:PIVOTER_RUN_REF"
  "ST:PIVOTER_RUN_ST"
  "RegionV2:PIVOTER_RUN_REGION_V2"
)

# ============ Helper: parse output ============
parse_result() {
  local output="$1"
  local algo="$2"

  # Total time
  total=$(echo "$output" | grep "NucleusCoreDecomposition took" | sed 's/.*took: //' | sed 's/ ms//' | head -1 || true)
  [ -z "$total" ] && total="N/A"

  # Peak memory
  mem=$(echo "$output" | grep "Final Memory" | grep -oE '[0-9]+' | tail -1 || true)
  [ -z "$mem" ] && mem="N/A"

  # Max core
  maxcore=$(echo "$output" | grep "Max core:" | sed 's/.*Max core: //' | head -1 || true)
  [ -z "$maxcore" ] && maxcore=$(echo "$output" | grep "^maxCoreVal=" | sed 's/.*=//' || true)
  [ -z "$maxcore" ] && maxcore="N/A"

  # Phase breakdown depends on algorithm
  enum_ms="N/A"; build_ms="N/A"; peel_ms="N/A"
  r_cliques="N/A"; r_tuples="N/A"; s_tuples="N/A"

  case "$algo" in
    RegionV2)
      enum_ms=$(echo "$output" | grep "MaxCliqEnum took" | sed 's/.*took: //' | sed 's/ ms//' || true)
      build_ms=$(echo "$output" | grep "Enumeration time:" | sed 's/.*: //' | sed 's/ ms//' || true)
      peel_ms=$(echo "$output" | grep "Peeling time:" | sed 's/.*: //' | sed 's/ ms//' || true)
      r_tuples=$(echo "$output" | grep "r-tuples:" | sed 's/.*r-tuples: //' | sed 's/,.*//' || true)
      s_tuples=$(echo "$output" | grep "s-tuples:" | sed 's/.*s-tuples: //' | head -1 || true)
      r_cliques=$(echo "$output" | grep "r-cliques:" | sed 's/.*r-cliques: //' | head -1 || true)
      ;;
    ST)
      enum_ms=$(echo "$output" | grep -E "SDCT_Fused took|SDCT_MaxClique took" | sed 's/.*took: //' | sed 's/ ms//' || true)
      build_ms=$(echo "$output" | grep "fused build+counting" | sed 's/.*took: //' | sed 's/ ms//' || true)
      peel_ms="N/A"
      ;;
    REF)
      enum_ms=$(echo "$output" | grep -E "SDCT_Fused took|SDCT_MaxClique took" | sed 's/.*took: //' | sed 's/ ms//' || true)
      build_ms=$(echo "$output" | grep "countingPerEdge" | sed 's/.*took: //' | sed 's/ ms//' || true)
      peel_ms="N/A"
      ;;
  esac

  [ -z "$enum_ms" ] && enum_ms="N/A"
  [ -z "$build_ms" ] && build_ms="N/A"
  [ -z "$peel_ms" ] && peel_ms="N/A"
  [ -z "$r_cliques" ] && r_cliques="N/A"
  [ -z "$r_tuples" ] && r_tuples="N/A"
  [ -z "$s_tuples" ] && s_tuples="N/A"
}

# ============ Verification helper ============
verify_result() {
  local graph=$1 r=$2 s=$3
  local st_out v2_out
  st_out=$(PIVOTER_RUN_ST=1 $BIN "graphs/${graph}.edges" "$r" "$s" 2>&1)
  v2_out=$(PIVOTER_RUN_REGION_V2=1 $BIN "graphs/${graph}.edges" "$r" "$s" 2>&1)

  st_cores=$(echo "$st_out" | grep "core=" | sed 's/^[[:space:]]*//' | sort -t= -k2 -n)
  v2_cores=$(echo "$v2_out" | grep "core=" | sed 's/^[[:space:]]*//' | sort -t= -k2 -n)

  if [ "$st_cores" = "$v2_cores" ]; then
    echo "EXACT"
  else
    echo "MISMATCH"
  fi
}

# ============ Main loop ============
run_count=0
total_runs=0
for ds in "${DATASETS[@]}"; do
  IFS=':' read -r graph omega <<< "$ds"
  for rs in $RS_COMBOS; do
    IFS=',' read -r r s <<< "$rs"
    [ "$s" -gt "$omega" ] && continue
    for alg in "${ALGOS[@]}"; do
      total_runs=$((total_runs + 1))
    done
  done
done

echo "Datasets: ${#DATASETS[@]}, (r,s) combos: $(echo $RS_COMBOS | wc -w | tr -d ' '), algorithms: ${#ALGOS[@]}"
echo "Total runs: $total_runs (timeout=${TIMEOUT}s)"
echo "Output: $OUTCSV"
echo ""

for ds in "${DATASETS[@]}"; do
  IFS=':' read -r graph omega <<< "$ds"

  if [ ! -f "graphs/${graph}.edges" ]; then
    echo "SKIP: graphs/${graph}.edges not found"
    continue
  fi

  echo "============================================================="
  echo "  Graph: $graph (ω=$omega)"
  echo "============================================================="

  for rs in $RS_COMBOS; do
    IFS=',' read -r r s <<< "$rs"
    [ "$s" -gt "$omega" ] && continue

    printf "\n  r=%s s=%s\n" "$r" "$s"
    printf "  %-12s %12s %10s %10s %10s %10s %10s\n" "Algorithm" "Total(ms)" "Enum(ms)" "Build(ms)" "Peel(ms)" "Mem(kB)" "MaxCore"
    printf "  %s\n" "$(printf -- '-%.0s' {1..82})"

    # Optional: verify first
    if [ $VERIFY -eq 1 ]; then
      vresult=$(verify_result "$graph" "$r" "$s" 2>/dev/null || echo "ERROR")
      printf "  Correctness: %s\n" "$vresult"
    fi

    for alg in "${ALGOS[@]}"; do
      IFS=':' read -r algo_name env_var <<< "$alg"
      run_count=$((run_count + 1))

      logfile="$LOGDIR/${graph}_r${r}_s${s}_${algo_name}.log"

      output=""
      exit_code=0
      output=$(env "${env_var}=1" $BIN "graphs/${graph}.edges" "$r" "$s" 2>&1) || exit_code=$?
      echo "$output" > "$logfile"

      status="OK"
      if [ $exit_code -eq 124 ]; then
        status="TIMEOUT"
        total="TIMEOUT"; enum_ms="N/A"; build_ms="N/A"; peel_ms="N/A"
        mem="N/A"; maxcore="N/A"; r_cliques="N/A"; r_tuples="N/A"; s_tuples="N/A"
      elif [ $exit_code -ne 0 ]; then
        status="ERROR($exit_code)"
        total="ERROR"; enum_ms="N/A"; build_ms="N/A"; peel_ms="N/A"
        mem="N/A"; maxcore="N/A"; r_cliques="N/A"; r_tuples="N/A"; s_tuples="N/A"
      else
        parse_result "$output" "$algo_name"
      fi

      echo "$graph,$r,$s,$algo_name,$total,$enum_ms,$build_ms,$peel_ms,$mem,$status,$maxcore,$r_cliques,$r_tuples,$s_tuples" >> "$OUTCSV"

      printf "  %-12s %12s %10s %10s %10s %10s %10s\n" \
        "$algo_name" "$total" "$enum_ms" "$build_ms" "$peel_ms" "$mem" "$maxcore"

    done
  done
  echo ""
done

# ============ Summary table ============
echo ""
echo "============================================================="
echo "  Summary: Speedup of RegionV2 vs REF and ST"
echo "============================================================="
echo ""
printf "%-18s %4s %4s | %12s %12s %12s | %8s %8s\n" \
  "Graph" "r" "s" "REF(ms)" "ST(ms)" "V2(ms)" "vs REF" "vs ST"
printf "%s\n" "$(printf -- '-%.0s' {1..96})"

# Parse CSV for summary
tail -n +2 "$OUTCSV" | awk -F',' '
{
  key = $1 "," $2 "," $3
  if ($4 == "REF") ref[key] = $5
  if ($4 == "ST") st[key] = $5
  if ($4 == "RegionV2") v2[key] = $5
  if (!(key in order)) { order[key] = NR; keys[NR] = key }
}
END {
  for (i = 1; i <= NR; i++) {
    if (!(i in keys)) continue
    key = keys[i]
    split(key, a, ",")
    r_t = (ref[key]+0 > 0) ? ref[key] : "N/A"
    s_t = (st[key]+0 > 0) ? st[key] : "N/A"
    v_t = (v2[key]+0 > 0) ? v2[key] : "N/A"
    vs_ref = "N/A"; vs_st = "N/A"
    if (v_t != "N/A" && v_t+0 > 0) {
      if (r_t != "N/A" && r_t+0 > 0) vs_ref = sprintf("%.1fx", r_t / v_t)
      if (s_t != "N/A" && s_t+0 > 0) vs_st = sprintf("%.1fx", s_t / v_t)
    }
    printf "%-18s %4s %4s | %12s %12s %12s | %8s %8s\n", a[1], a[2], a[3], r_t, s_t, v_t, vs_ref, vs_st
  }
}'

echo ""
echo "============================================================="
echo "  Memory Comparison"
echo "============================================================="
echo ""
printf "%-18s %4s %4s | %12s %12s %12s | %8s\n" \
  "Graph" "r" "s" "REF(kB)" "ST(kB)" "V2(kB)" "V2/ST"
printf "%s\n" "$(printf -- '-%.0s' {1..82})"

tail -n +2 "$OUTCSV" | awk -F',' '
{
  key = $1 "," $2 "," $3
  if ($4 == "REF") ref[key] = $9
  if ($4 == "ST") st[key] = $9
  if ($4 == "RegionV2") v2[key] = $9
  if (!(key in order)) { order[key] = NR; keys[NR] = key }
}
END {
  for (i = 1; i <= NR; i++) {
    if (!(i in keys)) continue
    key = keys[i]
    split(key, a, ",")
    r_m = (ref[key]+0 > 0) ? ref[key] : "N/A"
    s_m = (st[key]+0 > 0) ? st[key] : "N/A"
    v_m = (v2[key]+0 > 0) ? v2[key] : "N/A"
    ratio = "N/A"
    if (v_m+0 > 0 && s_m+0 > 0) ratio = sprintf("%.2fx", v_m / s_m)
    printf "%-18s %4s %4s | %12s %12s %12s | %8s\n", a[1], a[2], a[3], r_m, s_m, v_m, ratio
  }
}'

echo ""
echo "============================================================="
echo "  Region V2 Detail: Tuple Compression"
echo "============================================================="
echo ""
printf "%-18s %4s %4s | %12s %12s %12s | %10s\n" \
  "Graph" "r" "s" "r-cliques" "r-tuples" "s-tuples" "compress"
printf "%s\n" "$(printf -- '-%.0s' {1..80})"

tail -n +2 "$OUTCSV" | awk -F',' '$4=="RegionV2" {
  compress = "N/A"
  if ($13+0 > 0 && $12+0 > 0) compress = sprintf("%.0fx", $12 / $13)
  printf "%-18s %4s %4s | %12s %12s %12s | %10s\n", $1, $2, $3, $12, $13, $14, compress
}'

echo ""
echo "============================================================="
echo "  DONE. Results: $OUTCSV"
echo "  Logs: $LOGDIR/"
echo "============================================================="
