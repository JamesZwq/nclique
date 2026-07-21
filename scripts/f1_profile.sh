exec >/data/wenqianz/f1_profile.out 2>&1
set -u
echo "=== F1 DECISION EXPERIMENT START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter
GD=/data/wenqianz
cd $REPO
# server tree was verified identical to origin/main except for the F1 block, so reset is safe
git fetch -q origin && git reset -q --hard origin/main
echo "HEAD now: $(git log --oneline -1)"
grep -c F1_PROFILE region_native/region_native_sct_peel.cpp | sed 's/^/F1_PROFILE occurrences in source: /'
cd $REPO/region_native
echo "--- rebuilding rn_diag ---"
g++ -O3 -march=native -std=c++17 -fopenmp -I. -I../src/NucleusDecomposition \
    -o rn_diag region_native_sct_peel.cpp 2>&1 | head -20
echo "build_exit=$?"
export SCT_MAX_INC=200000000
RN=$REPO/region_native/rn_diag

for G in web-Google ca-AstroPh; do
  GR=$GD/$G.edges; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  echo ""; echo "########## $G (plane r=3..5, s<=8, 1 thread) ##########"
  timeout 3000 /usr/bin/time -v env OMP_NUM_THREADS=1 SCT_RSWEEP=1 \
      SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 F1_PROFILE=1 \
      $RN $GR 3 4 >/dev/null 2>/tmp/f1_$G
  echo "exit=$?"
  grep -E '^\[f1\]' /tmp/f1_$G | sed "s/^/RAW $G /"
  grep -E 'Elapsed \(wall|Maximum resident' /tmp/f1_$G | sed 's/^/  /'
done
echo ""; echo "=== F1_DONE $(date) ==="
