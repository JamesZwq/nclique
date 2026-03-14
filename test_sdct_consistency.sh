#!/bin/bash

GRAPH="./graphs/com-dblp.edges"
BIN="./build-opt/bin/degeneracy_cliques"

echo "Testing SDCT vs SDCT_Par consistency..."
echo "Graph: $GRAPH"
echo ""

# Test with small graph first
export PIVOTER_RUN_REF=True
export OMP_NUM_THREADS=4

echo "Running with 4 threads..."
$BIN $GRAPH 3 4 degen 2>&1 | tee /tmp/sdct_test.log

echo ""
echo "Checking results..."
grep "Core Value Distribution:" -A 50 /tmp/sdct_test.log | head -20
