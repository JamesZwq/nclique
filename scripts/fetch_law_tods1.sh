#!/bin/bash
# Fetch LAW billion-scale web graphs on tods1 (NOT tods2: it is timing FEM scouts; tods1 is the
# high-load workhorse where contention does not matter). WebGraph BVGraph -> undirected .edges.
# Pipeline per graph: .graph+.properties -> regenerate offsets -> ASCII adjacency dump ->
# symmetrize (min,max) -> sort -u -> "n m" header + pairs.
OUT=/data/wenqianz/fetch_law.out
exec >"$OUT" 2>&1
set -u
D=/data/wenqianz/law; mkdir -p "$D"; cd "$D" || exit 1
echo "=== FETCH LAW START $(date) ==="

# --- webgraph tooling (site bundle includes all deps) ---
if [ ! -d deps ] || ! ls webgraph-*.jar >/dev/null 2>&1; then
  ok=0
  for V in 3.6.11 3.6.10 3.6.9 3.6.8; do
    wget -q "https://webgraph.di.unimi.it/webgraph-$V-bin.tar.gz" -O wg.tar.gz || continue
    tar xzf wg.tar.gz && cp webgraph-$V/webgraph-$V.jar . && ok=1 && echo "webgraph $V" && break
  done
  [ $ok -eq 1 ] || { echo "FATAL: webgraph jar download failed (all versions tried)"; exit 1; }
  wget -q "https://webgraph.di.unimi.it/webgraph-deps.tar.gz" -O deps.tar.gz \
    || { echo "FATAL: webgraph-deps download failed"; exit 1; }
  mkdir -p deps && tar xzf deps.tar.gz -C deps
  rm -f wg.tar.gz deps.tar.gz
fi
CP="$(ls webgraph-*.jar | head -1):$(ls deps/*.jar | tr '\n' ':')"
echo "classpath jars: $(ls deps/*.jar | wc -l) + main"

for G in it-2004 uk-2005; do
  [ -f "/data/wenqianz/$G.full.edges" ] && { echo "$G.full.edges exists, skip"; continue; }
  echo "--- $G ---"
  for ext in properties graph; do
    [ -f "$G.$ext" ] && continue
    wget -q "http://data.law.di.unimi.it/webdata/$G/$G.$ext" -O "$G.$ext" \
      || { echo "  DOWNLOAD FAILED $G.$ext"; rm -f "$G.$ext"; continue 2; }
  done
  echo "  graph=$(ls -la $G.graph | awk '{print $5}')B"
  # offsets (-o -O -L regenerates offsets + load test)
  java -Xmx100g -cp "$CP" it.unimi.dsi.webgraph.BVGraph -o -O -L "$G" \
    || { echo "  OFFSETS FAILED $G"; continue; }
  # ASCII dump: ${G}-ascii.graph-txt = line1 n, then successor list per node
  java -Xmx100g -cp "$CP" it.unimi.dsi.webgraph.ASCIIGraph "$G" "$G-ascii" \
    || { echo "  ASCII DUMP FAILED $G"; continue; }
  # symmetrize + dedupe (polite parallelism: the box is shared)
  awk 'NR==1 { n=$1; next }
       { u=NR-2; for (i=1;i<=NF;i++) { v=$i; if (v==u) continue;
           if (u<v) print u, v; else print v, u } }' "$G-ascii.graph-txt" \
    | sort -u -S 40G --parallel=8 -T "$D" > "$G.pairs"
  M=$(wc -l < "$G.pairs")
  N=$(head -1 "$G-ascii.graph-txt")
  { echo "$N $M"; cat "$G.pairs"; } > "/data/wenqianz/$G.full.edges"
  rm -f "$G.pairs" "$G-ascii.graph-txt"
  echo "  DONE: /data/wenqianz/$G.full.edges  header=$N $M  size=$(ls -la /data/wenqianz/$G.full.edges | awk '{print $5}')B"
  df -h /data | tail -1
done
echo "=== FETCH LAW DONE $(date) ==="
