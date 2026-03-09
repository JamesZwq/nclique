#!/bin/bash

echo "=== Building degeneracy_cliques ==="
cd /home/wenqianz/nclique/build
make -j8 degeneracy_cliques

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo ""
echo "=== Testing com-dblp (smaller dataset) ==="
./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep -E '(Multi-threading|Threads:|speedup|took:)'

echo ""
echo "=== Testing web-Google ==="
./bin/degeneracy_cliques /data/wenqianz/web-Google.edges 3 4 2>&1 | grep -E '(Multi-threading|Threads:|speedup|took:)'

echo ""
echo "=== Test completed ==="
