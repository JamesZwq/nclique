exec >/data/wenqianz/step2_headroom.out 2>&1
set -u
echo "=== STEP2 HEADROOM PROBE START $(date) ==="
echo "Measures witnessProbes vs newWitnesses. probes/new = the redundancy factor that a"
echo "grouped walk could remove. If ~1.0 there is NOTHING to gain and step 2 is dead."
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
cd $REPO && git fetch -q origin && git reset -q --hard origin/main
echo "HEAD: $(git log --oneline -1)"
cd $REPO/region_native
[ -x rn_slp ] || g++ -O3 -march=native -std=c++17 -fopenmp -I. -I../src/NucleusDecomposition -o rn_slp region_native_sct_peel.cpp
echo "binary=$(ls -la rn_slp | awk '{print $5}')"
export SCT_MAX_INC=200000000
for G in ca-AstroPh web-Google; do
  GR=$GD/$G.edges; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  echo ""; echo "########## $G plane r=3..5 s<=8 ##########"
  timeout 3000 env OMP_NUM_THREADS=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
      PIVOTER_PEEL_PROFILE=1 ./rn_slp $GR 3 4 >/dev/null 2>/tmp/h_$G; echo "exit=$?"
  grep -E 'witness-delta' /tmp/h_$G | sed 's/^/  /'
done
echo ""; echo "=== STEP2_HEADROOM_DONE $(date) ==="
