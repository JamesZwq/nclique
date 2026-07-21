cd /Users/zhangwenqian/UNSW/pivoter/region_native
PASS=0; FAIL=0
run(){ # run <bin> <graph> <rmin> <rmax> <smax> <out>
  OMP_NUM_THREADS=1 SCT_RSWEEP=1 SCT_RMIN=$3 SCT_RMAX=$4 SCT_SMAX=$5 \
    "$1" "$2" "$3" $(( $3 + 1 )) 2>/dev/null \
    | grep -E '^core=|^\[sct-peel\]|^\[nsi-plane-cell\]' \
    | sed -E 's/replay=[0-9.]+s//' > "$6"
}
for G in .multir_gate.edges .multir_residue_gate.edges ../raw_datasets/dblp-sigmod.edges ../nuclearSigmod/toy.edges; do
  [ -f "$G" ] || continue
  for cfg in "3 5 8" "3 4 6" "2 4 7" "3 3 8"; do
    rmin=$(echo $cfg | cut -d' ' -f1); rmax=$(echo $cfg | cut -d' ' -f2); smax=$(echo $cfg | cut -d' ' -f3)
    run /tmp/rn_old  "$G" $rmin $rmax $smax /tmp/o.txt
    run /tmp/rn_slp1 "$G" $rmin $rmax $smax /tmp/n.txt
    n=$(wc -l < /tmp/o.txt)
    if [ "$n" -gt 0 ] && diff -q /tmp/o.txt /tmp/n.txt >/dev/null 2>&1; then
      echo "  PASS $(basename $G) r=$rmin..$rmax s<=$smax ($n lines)"; PASS=$((PASS+1))
    elif [ "$n" -eq 0 ]; then
      echo "  SKIP $(basename $G) r=$rmin..$rmax s<=$smax (no output - invalid config)"
    else
      echo "  *** FAIL $(basename $G) r=$rmin..$rmax s<=$smax"; diff /tmp/o.txt /tmp/n.txt | head -8; FAIL=$((FAIL+1))
    fi
  done
done
echo "=== BIT-EXACT GATE: PASS=$PASS FAIL=$FAIL ==="
