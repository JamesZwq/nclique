exec >/data/wenqianz/slp_measure.out 2>&1
set -u
echo "=== SLP STEP1 + F2 MEASUREMENT START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
cd $REPO && git fetch -q origin && git reset -q --hard origin/main
echo "HEAD: $(git log --oneline -1)"
echo "guard: liveResidueLeaf occurrences = $(grep -c liveResidueLeaf region_native/region_native_sct_peel.cpp)"
echo "guard: F2_PROFILE occurrences      = $(grep -c F2_PROFILE region_native/region_native_sct_peel.cpp)"
cd $REPO/region_native
g++ -O3 -march=native -std=c++17 -fopenmp -I. -I../src/NucleusDecomposition -o rn_slp region_native_sct_peel.cpp 2>&1 | head -10
echo "build_new=$?"
git stash -q 2>/dev/null
git show HEAD~1:region_native/region_native_sct_peel.cpp > /tmp/old_sct.cpp
g++ -O3 -march=native -std=c++17 -fopenmp -I. -I../src/NucleusDecomposition -o rn_base /tmp/old_sct.cpp 2>&1 | head -10
echo "build_base=$?"
export SCT_MAX_INC=200000000

for G in ca-AstroPh web-Google; do
  GR=$GD/$G.edges; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  echo ""; echo "########## $G plane r=3..5 s<=8 ##########"
  echo "--- BASELINE (pre-§201) ---"
  timeout 3000 /usr/bin/time -v env OMP_NUM_THREADS=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
      ./rn_base $GR 3 4 >/tmp/base_$G.out 2>/tmp/base_$G.err; echo "exit=$?"
  grep -E 'Elapsed \(wall|Maximum resident' /tmp/base_$G.err | sed 's/^/  /'
  echo "--- SLP STEP1 (with F1+F2 profiling) ---"
  timeout 3000 /usr/bin/time -v env OMP_NUM_THREADS=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
      F1_PROFILE=1 F2_PROFILE=1 ./rn_slp $GR 3 4 >/tmp/slp_$G.out 2>/tmp/slp_$G.err; echo "exit=$?"
  grep -E 'Elapsed \(wall|Maximum resident' /tmp/slp_$G.err | sed 's/^/  /'
  echo "  --- BIT-EXACTNESS vs baseline ---"
  grep -E '^core=|^\[sct-peel\]' /tmp/base_$G.out > /tmp/b.f; grep -E '^core=|^\[sct-peel\]' /tmp/slp_$G.out > /tmp/s.f
  if diff -q /tmp/b.f /tmp/s.f >/dev/null 2>&1; then echo "  BIT-EXACT: PASS ($(wc -l < /tmp/b.f) lines)"; else echo "  BIT-EXACT: *** FAIL ***"; diff /tmp/b.f /tmp/s.f | head -10; fi
  echo "  --- F1 (certified-replay share) + SLP skips ---"
  grep -E '^\[f1\]|^\[slp\]' /tmp/slp_$G.err | sed 's/^/  /'
  echo "  --- F2 (T3+ would-add) ---"
  grep -E '^\[f2\]' /tmp/slp_$G.err | sed 's/^/  /'
done
echo ""; echo "=== SLP_MEASURE_DONE $(date) ==="
