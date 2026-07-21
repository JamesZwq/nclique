exec >/data/wenqianz/astro_recompare.out 2>&1
set -u
echo "=== ca-AstroPh OURS(post-§201) vs CND RECOMPARE START $(date) ==="
echo "Does the SLP step-1 2.01x flip the cells where CND used to beat us?"
REPO=/home/wenqianz/UNSW/pivoter; GD=/data/wenqianz
cd $REPO && git fetch -q origin && git reset -q --hard origin/main
echo "HEAD: $(git log --oneline -1)"
cp $REPO/pivoter/nCr.txt $REPO/src/nCr.txt 2>/dev/null
cd $REPO/region_native
[ -x rn_slp ] || g++ -O3 -march=native -std=c++17 -fopenmp -I. -I../src/NucleusDecomposition -o rn_slp region_native_sct_peel.cpp
RN=$REPO/region_native/rn_slp; CND=$REPO/build/bin/degeneracy_cliques
export SCT_MAX_INC=200000000
GR=$GD/ca-AstroPh.edges
tsec(){ grep -i 'Elapsed (wall' "$1" 2>/dev/null|awk '{v=$NF;n=split(v,a,":");if(n==3)print a[1]*3600+a[2]*60+a[3];else print a[1]*60+a[2]}'; }
rss(){ awk -F': ' '/Maximum resident/{printf "%.1f",$2/1048576}' "$1" 2>/dev/null; }

echo ""; echo "########## A. WHOLE-PLANE (ours, ONE build, r=3..5 s<=8, 1 thread) ##########"
/usr/bin/time -v env OMP_NUM_THREADS=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 \
    $RN $GR 3 4 >/dev/null 2>/tmp/ac_plane
echo "  OUR WHOLE PLANE: $(tsec /tmp/ac_plane)s  $(rss /tmp/ac_plane)GB"

echo ""; echo "########## B. PER-CELL: ours (fixed-r sweep) vs CND (96-core) ##########"
CTO=1800
for r in 3 4 5; do for s in $(seq $((r+1)) 8); do
  # CND single cell
  ( cd $REPO && timeout $CTO /usr/bin/time -v $CND $GR $r $s >/dev/null 2>/tmp/ac_c_${r}_${s} ); cex=$?
  ct=$(tsec /tmp/ac_c_${r}_${s}); cm=$(rss /tmp/ac_c_${r}_${s})
  [ $cex -ne 0 ] && ct="FAIL(exit=$cex)"
  echo "  ($r,$s)  CND(96c): ${ct}s ${cm}GB"
done; done
echo ""; echo "=== ASTRO_RECOMPARE_DONE $(date) ==="
