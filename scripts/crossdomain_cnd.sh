exec >/data/wenqianz/crossdomain_cnd.out 2>&1
set -u
echo "=== CROSSDOMAIN PLANE_VS_CND START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter
RN=$REPO/region_native/rn_diag; CND=$REPO/build/bin/degeneracy_cliques; GD=/data/wenqianz
cp $REPO/pivoter/nCr.txt $REPO/src/nCr.txt 2>/dev/null
export SCT_MAX_INC=200000000
PTO=2400; CTO=1200    # plane / per-CND-cell timeouts (tightened to bound total runtime)
tsec(){ grep -i 'Elapsed (wall' "$1" 2>/dev/null|awk '{v=$NF;n=split(v,a,":");if(n==3)print a[1]*3600+a[2]*60+a[3];else print a[1]*60+a[2]}'; }
rss(){ awk -F': ' '/Maximum resident/{printf "%.1f",$2/1048576}' "$1" 2>/dev/null; }
cd $REPO/region_native
[ -x rn_diag ] || g++ -O3 -march=native -std=c++17 -fopenmp -I. -I../src/NucleusDecomposition -o rn_diag region_native_sct_peel.cpp
# domains (non-collaboration): co-purchase, web (small->large), email. cleanest-completing FIRST.
for G in com-amazon.ungraph web-Stanford web-Google email-Eu-core; do
  GR=$GD/$G.edges; [ -f "$GR" ] || { echo "$G MISSING"; continue; }
  echo "===== $G ====="
  timeout $PTO /usr/bin/time -v env OMP_NUM_THREADS=1 SCT_RSWEEP=1 SCT_RMIN=3 SCT_RMAX=5 SCT_SMAX=8 $RN $GR 3 4 >/dev/null 2>/tmp/pl_$G; pex=$?
  pt=$(tsec /tmp/pl_$G); pm=$(rss /tmp/pl_$G); [ $pex -ne 0 ] && pt="TO/OOM(exit=$pex)"
  echo "  OUR PLANE (1-thread, ONE build, r=3..5 s<=8): ${pt}s ${pm}GB"
  cs=0; nc=0; fc=0; oom=""
  for r in 3 4 5; do for s in $(seq $((r+1)) 8); do
    ( cd $REPO && timeout $CTO /usr/bin/time -v $CND $GR $r $s >/dev/null 2>/tmp/c_${G}_${r}_${s} ); ex=$?
    t=$(tsec /tmp/c_${G}_${r}_${s}); m=$(rss /tmp/c_${G}_${r}_${s})
    if [ $ex -eq 0 ] && [ -n "$t" ]; then cs=$(awk -v a=$cs -v b=$t 'BEGIN{print a+b}'); nc=$((nc+1)); echo "    CND ($r,$s): ${t}s ${m}GB"
    else fc=$((fc+1)); [ $ex -eq 124 ] && oom="$oom ($r,$s)=TO" || oom="$oom ($r,$s)=exit$ex(rss${m}GB)"; echo "    CND ($r,$s): FAIL exit=$ex rss=${m}GB"; fi
  done; done
  echo "  CND(96-core) SUM over $nc cells ($fc failed):$oom  => ${cs}s"
  [ -n "$pt" ] && echo "  >>> $G: PLANE ${pt}s vs CND ${cs}s/$nc-cells => $(awk -v c=$cs -v p=$pt 'BEGIN{if(p+0>0)printf "%.1f",c/p;else print "NA"}')x (+ $fc CND cells infeasible)"
  echo "===== $G DONE $(date) ====="
done
echo "=== CROSSDOMAIN_DONE $(date) ==="
