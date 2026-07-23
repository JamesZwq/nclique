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

# --- webgraph tooling via Maven Central (the site tarball path 404s; maven URLs are stable) ---
M=https://repo1.maven.org/maven2
JARS="
it/unimi/dsi/webgraph/3.6.10/webgraph-3.6.10.jar
it/unimi/dsi/dsiutils/2.7.3/dsiutils-2.7.3.jar
it/unimi/dsi/fastutil/8.5.12/fastutil-8.5.12.jar
it/unimi/dsi/sux4j/5.2.3/sux4j-5.2.3.jar
com/martiansoftware/jsap/2.1/jsap-2.1.jar
org/slf4j/slf4j-api/1.7.36/slf4j-api-1.7.36.jar
org/slf4j/slf4j-simple/1.7.36/slf4j-simple-1.7.36.jar
org/apache/commons/commons-lang3/3.12.0/commons-lang3-3.12.0.jar
commons-io/commons-io/2.11.0/commons-io-2.11.0.jar
org/apache/commons/commons-configuration2/2.8.0/commons-configuration2-2.8.0.jar
commons-beanutils/commons-beanutils/1.9.4/commons-beanutils-1.9.4.jar
org/apache/commons/commons-collections4/4.4/commons-collections4-4.4.jar
commons-logging/commons-logging/1.2/commons-logging-1.2.jar
org/apache/commons/commons-text/1.10.0/commons-text-1.10.0.jar
com/google/guava/guava/31.1-jre/guava-31.1-jre.jar
"
mkdir -p deps
for j in $JARS; do
  f=deps/$(basename "$j")
  [ -f "$f" ] && continue
  wget -q "$M/$j" -O "$f" || { echo "FATAL: maven download failed: $j"; exit 1; }
done
CP="$(ls deps/*.jar | tr '\n' ':')"
echo "classpath jars: $(ls deps/*.jar | wc -l)"

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
