#!/bin/bash
# Fetch FEM/mesh matrices from SuiteSparse, convert MatrixMarket -> .edges, validate each.
# Old roster graphs are re-fetched because their files survive on NEITHER server (only logs), and
# their old "CND infeasible" claims are §217-quarantined anyway (parallel-CND era) -- they must be
# re-earned serially. New big ones extend the roster toward the acceptance standard.
# .edges format: first line "n m", then one "u v" pair per line, 0-based, each undirected edge once.
OUT=/data/wenqianz/fetch_fem.out
exec >"$OUT" 2>&1
set -u
cd /data/wenqianz || exit 1
echo "=== FETCH FEM START $(date) ==="

fetch() { # $1=group $2=name
  local g=$1 n=$2
  [ -f "$n.edges" ] && { echo "$n.edges exists, skip"; return 0; }
  echo "--- $g/$n ---"
  wget -q "https://suitesparse-collection-website.herokuapp.com/MM/$g/$n.tar.gz" -O "$n.tar.gz" \
    || { echo "  DOWNLOAD FAILED $n"; rm -f "$n.tar.gz"; return 1; }
  tar xzf "$n.tar.gz" || { echo "  UNTAR FAILED $n"; return 1; }
  # convert: symmetric MatrixMarket, 1-based, lower-triangle; emit 0-based unique undirected edges
  awk -v out="$n.edges.tmp" '
    /^%/ { next }
    !dims { n=$1; dims=1; next }
    { i=$1-1; j=$2-1; if (i==j) next;
      if (i<j) { t=i; i=j; j=t }        # normalize i>j
      print i, j > out; m++ }
    END { printf "%d %d\n", n, m > (out ".hdr") }
  ' "$n/$n.mtx" || { echo "  CONVERT FAILED $n"; return 1; }
  cat "$n.edges.tmp.hdr" "$n.edges.tmp" > "$n.edges"
  rm -rf "$n.tar.gz" "$n" "$n.edges.tmp" "$n.edges.tmp.hdr"
  echo "  $n.edges: $(head -1 "$n.edges")  ($(ls -la "$n.edges" | awk '{print $5}') bytes)"
}

# validation first (small), then the roster, then the big ones
fetch Nasa      nasasrb    && VAL=nasasrb
fetch Boeing    pwtk
fetch Chen      pkustk11
fetch Chen      pkustk13
fetch GHS_psdef ldoor
fetch Janna     Flan_1565
fetch Janna     Queen_4147

# pipeline validation: the engine must load nasasrb and survive a (3,4) front-end
if [ -n "${VAL:-}" ]; then
  echo "--- pipeline validation on $VAL ---"
  cd /home/wenqianz/UNSW/pivoter/region_native \
    && g++ -O3 -march=native -std=c++17 -I. -I../src/NucleusDecomposition -o /tmp/rn_val region_native_sct_peel.cpp \
    && SCT_W_ONLY=1 SCT_MAX_INC=500000000 /tmp/rn_val /data/wenqianz/$VAL.edges 3 4 2>&1 | grep -E '\[rn\]|\[W ' | head -6
fi
df -h /data | tail -1
echo "=== FETCH FEM DONE $(date) ==="
