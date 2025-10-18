#!/usr/bin/env bash
#  degeneracy_cliques， /bin/time -v，：experimentdata
# ：chmod +x run_all.sh && ./run_all.sh
set -u -o pipefail

TIME_BIN="/bin/time"   #  /bin/time -v
if [[ ! -x "$TIME_BIN" ]]; then
  echo "[FATAL] /bin/time " >&2
  exit 1
fi

BIN="./build/bin/degeneracy_cliques"
if [[ ! -x "$BIN" ]]; then
  echo "[FATAL] ：$BIN" >&2
  exit 1
fi

# （）
OUTPUT_FILE="experimentdata"

# 
DATASETS=(
  "/data/_/web-it-2004.edges"
  "/data/_/com-dblp.edges"
  "/data/_/web-Google.edges"
  "/data/_/soc-pokec-relationships.edges"
  "/data/_/web-Stanford.edges"
  "/data/_/com-youtube.edges"
)

# s  r ：
# s=1 -> r=3..30
# s=2 -> r=4..30
# s=3 -> r=4..30
S_VALUES=(3 4)

r_range_for_s () {
  local s="$1"
  case "$s" in
    1) echo "3 30" ;;
    2|3) echo "4 30" ;;
    *) echo "5 30" ;; # ，
  esac
}

TMPDIR="$(mktemp -d)"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT INT TERM

echo "===== RUN START $(date '+%F %T') =====" >> "$OUTPUT_FILE"

for ds in "${DATASETS[@]}"; do
  base="$(basename "$ds")"
  if [[ ! -f "$ds" ]]; then
    {
      echo
      echo "----- $(date '+%F %T') DATASET NOT FOUND -----"
      echo "dataset: $ds"
    } >> "$OUTPUT_FILE"
    echo "[WARN] dataset not found: $ds" >&2
    continue
  fi

  for s in "${S_VALUES[@]}"; do
    read -r RSTART REND < <(r_range_for_s "$s")
    for (( r=RSTART; r<=REND; ++r )); do
      ts="$(date '+%F %T')"
      echo "[RUN] $base s=$s r=$r @ $ts"

      {
        echo
        echo "================== RUN =================="
        echo "timestamp : $ts"
        echo "dataset   : $ds"
        echo "s         : $s"
        echo "r         : $r"
        echo "binary    : $BIN"
        echo "time_bin  : $TIME_BIN -v"
        echo "cmd       : $BIN $ds $s $r"
        echo "---------------- PROGRAM OUTPUT ----------------"
      } >> "$OUTPUT_FILE"

      tlog="$TMPDIR/time.$$.$RANDOM.log"
      #  stdout/stderr  OUTPUT_FILE；/bin/time -v ，
      "$TIME_BIN" -v -o "$tlog" "$BIN" "$ds" "$s" "$r" >> "$OUTPUT_FILE" 2>&1
      status=$?

      {
        echo "---------------- /bin/time -v ----------------"
        cat "$tlog"
        echo "----------------   STATUS   ------------------"
        echo "exit_status: $status"
        echo "================ END OF RUN =================="
      } >> "$OUTPUT_FILE"

      rm -f "$tlog"
    done
  done
done

echo "===== RUN END $(date '+%F %T') =====" >> "$OUTPUT_FILE"
echo "：$OUTPUT_FILE"