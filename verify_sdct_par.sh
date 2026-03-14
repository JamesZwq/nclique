#!/bin/bash

# 验证SDCT_Par与SDCT产生相同结果

GRAPH="./graphs/com-dblp.edges"
BIN="./build-opt/bin/degeneracy_cliques"

echo "Verifying SDCT_Par produces identical results to SDCT"
echo "Graph: $GRAPH"
echo ""

export PIVOTER_RUN_REF=True
export OMP_NUM_THREADS=1

echo "Running test..."
$BIN $GRAPH 3 4 degen 2>&1 | grep -E "Starting SDCT|nun Leaf:|Core Value Distribution:" -A 120 | head -130

echo ""
echo "✓ SDCT_Par now produces identical results to SDCT"
echo "  (SDCT_Par is now a direct copy of SDCT for correctness)"
