#!/bin/bash
cd /home/wenqianz/nclique

echo '========================================'
echo 'Baseline Performance Test'
echo 'Date: '$(date)
echo '========================================'
echo ''

# 测试单线程 baseline
echo 'Testing SINGLE THREAD (baseline)...'
export OMP_NUM_THREADS=1
for run in 1 2 3; do
  echo "  Run $run:"
  timeout 300 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep -E '(Multi-threading|took:)'
  sleep 2
done

echo ''

# 测试 16 线程
echo 'Testing 16 THREADS...'
export OMP_NUM_THREADS=16
for run in 1 2 3; do
  echo "  Run $run:"
  timeout 300 ./build/bin/degeneracy_cliques /data/wenqianz/com-dblp.edges 3 4 2>&1 | grep -E '(Multi-threading|took:)'
  sleep 2
done

echo ''
echo 'Test completed!'
