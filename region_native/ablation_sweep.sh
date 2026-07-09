#!/usr/bin/env bash
# §162 THREE-ALGORITHM ABLATION for the NSI sweep (base + 2 optimizations).
# Runs, per graph, three configs of region_native_sct_peel at UNIFORM r=4, s=5..8:
#   base : SpecND-Base        SCT_NO_RMERGE=1 SCT_SWEEP_NOCERT=1   (shared CPI, per-cell FULL peel)
#   cf   : + closed-form regs SCT_SWEEP_NOCERT=1                    (mergeable ON, certify OFF)
#   full : + certification    (default)                            (both optimizations ON = SpecND)
# CLEAN TIMING: single-thread (OMP=1), SERIAL within a machine. Run two machines in PARALLEL,
# each machine its own disjoint graph list. NEVER run two configs concurrently on one box.
# Metrics: wall (elapsed) + peak RSS from /usr/bin/time -v; certified% + mergeable from stdout.
#
# Usage: ./ablation_sweep.sh <machine_label> <bin> <graph_dir> <out_tsv> <graph1> [graph2 ...]
#   graphN = basename WITHOUT .edges (e.g. ca-AstroPh)
set -u
LABEL="$1"; BIN="$2"; GDIR="$3"; OUT="$4"; shift 4
GRAPHS=("$@")
R=4; S0=5; SMAX=8
TIMEOUT="${ABL_TIMEOUT:-1800}"          # per-config wall cap (s); a config that can't finish is infeasible
LOGDIR="$(dirname "$OUT")/ablog_${LABEL}"
mkdir -p "$LOGDIR"
TIMEBIN="$(command -v /usr/bin/time || echo time)"

if [ ! -f "$OUT" ]; then
  printf "machine\tgraph\tconfig\tr\ts0\tsmax\tstatus\twall_s\tpeak_rss_kb\tmergeable_rc\tcert_pct_last\texit\n" > "$OUT"
fi

run_one () {
  local g="$1" cfg="$2" env_toggles="$3"
  local gf="$GDIR/$g.edges"
  local stdout="$LOGDIR/${g}__${cfg}.out" tlog="$LOGDIR/${g}__${cfg}.time"
  if [ ! -f "$gf" ]; then
    printf "%s\t%s\t%s\t%d\t%d\t%d\tNOGRAPH\t\t\t\t\t\n" "$LABEL" "$g" "$cfg" "$R" "$S0" "$SMAX" >> "$OUT"
    echo "[skip] $g $cfg: no graph file $gf"; return
  fi
  echo "[run ] $LABEL $g $cfg  (timeout ${TIMEOUT}s)"
  # /usr/bin/time WRAPS timeout so RSS+elapsed are reported even on a kill.
  env SCT_SWEEP="$SMAX" OMP_NUM_THREADS=1 $env_toggles \
    "$TIMEBIN" -v timeout "$TIMEOUT" "$BIN" "$gf" "$R" "$S0" \
    > "$stdout" 2> "$tlog"
  local ec=$?
  local wall rss status
  wall=$(grep -a "Elapsed (wall clock)" "$tlog" | grep -aoE "[0-9:.]+$" | head -1 | awk -F: '{if(NF==3)print $1*3600+$2*60+$3; else if(NF==2)print $1*60+$2; else print $1}')
  rss=$(grep -a "Maximum resident set size" "$tlog" | grep -aoE "[0-9]+" | head -1)
  local mrg cert
  mrg=$(grep -a "r-mergeable:" "$stdout" | grep -aoE "\([0-9]+ r-cliques\)" | grep -aoE "[0-9]+" | head -1)
  cert=$(grep -a "chain-certified=" "$stdout" | tail -1 | grep -aoE "\([0-9.]+% of r-cliques\)" | grep -aoE "[0-9.]+" | head -1)
  if [ "$ec" -eq 124 ]; then status="TIMEOUT"; elif [ "$ec" -eq 0 ]; then status="OK"; else status="ERR$ec"; fi
  printf "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\n" \
    "$LABEL" "$g" "$cfg" "$R" "$S0" "$SMAX" "$status" "${wall:-}" "${rss:-}" "${mrg:-}" "${cert:-}" "$ec" >> "$OUT"
  echo "       -> $status  wall=${wall:-?}s  rss=${rss:-?}kB  mergeable=${mrg:-?}  cert%=${cert:-?}"
}

echo "=== ablation start $LABEL  $(date)  graphs: ${GRAPHS[*]} ==="
for g in "${GRAPHS[@]}"; do
  # fast config first (full SpecND) so we always capture the headline; then the two slow ones.
  run_one "$g" full ""
  run_one "$g" cf   "SCT_SWEEP_NOCERT=1"
  run_one "$g" base "SCT_NO_RMERGE=1 SCT_SWEEP_NOCERT=1"
done
echo "=== ablation DONE $LABEL  $(date) ==="
echo "ABLATION_RC=0 $LABEL"
