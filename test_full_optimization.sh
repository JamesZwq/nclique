#!/bin/bash

echo "=== Building with full optimizations ==="
cd /home/wenqianz/nclique/build
make -j8 degeneracy_cliques

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo ""
echo "=== Testing com-dblp ==="
echo "Running with optimized SDCT_Par + optimized Nucleus..."
./bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | tee com-dblp-result.txt

echo ""
echo "=== Testing web-Google ==="
echo "Running with optimized SDCT_Par + optimized Nucleus..."
./bin/degeneracy_cliques /data/wenqianz/web-Google.edges 3 4 2>&1 | tee web-google-result.txt

echo ""
echo "=== Extracting performance metrics ==="
echo ""
echo "com-dblp results:"
grep -E "(Tree Build took:|Multi-threading|speedup)" com-dblp-result.txt

echo ""
echo "web-Google results:"
grep -E "(Tree Build took:|Multi-threading|speedup)" web-google-result.txt

echo ""
echo "=== Test completed ==="
