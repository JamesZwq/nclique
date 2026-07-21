exec >/data/wenqianz/email_diag.out 2>&1
set -u
echo "=== EMAIL-vs-AMAZON DIAGNOSTIC START $(date) ==="
REPO=/home/wenqianz/UNSW/pivoter
RN=$REPO/region_native/rn_diag; GD=/data/wenqianz
export SCT_MAX_INC=200000000
cd $REPO/region_native

tsec(){ grep -i 'Elapsed (wall' "$1" 2>/dev/null|awk '{v=$NF;n=split(v,a,":");if(n==3)print a[1]*3600+a[2]*60+a[3];else print a[1]*60+a[2]}'; }

# --- fixed-r SWEEP mode: this path REACHES the compression counters (plane mode returns at 1423) ---
sweep(){ # sweep <tag> <graph> <r> <smax>
  tag=$1; g=$2; r=$3; sm=$4; L=/tmp/ed_$tag
  timeout 2400 /usr/bin/time -v env OMP_NUM_THREADS=1 SCT_SWEEP=$sm \
    REGCLS_DBG=1 COMP_DBG=1 CLS_LEAF_DBG=1 PE_DBG=1 MAPS_DBG=1 MEM_DBG=1 RM_DBG=1 \
    $RN $GD/$g.edges $r $((r+1)) >/dev/null 2>$L; ex=$?
  echo "  [$tag] $g fixed-r=$r smax=$sm -> $(tsec $L)s (exit=$ex)"
  grep -E '^\[(comp|cls-leaf|regcls|pe-dbg|maps-dbg)\]|PATTERN-EXPLOSION' $L 2>/dev/null | sed 's/^/      /' | head -12
  grep -E 'chain-certified' $L 2>/dev/null | sed 's/^/      /' | head -8
}

# --- plane mode: timing decomposition only ---
plane(){ # plane <tag> <graph> <rmin> <rmax> <smax>
  tag=$1; g=$2; rmin=$3; rmax=$4; sm=$5; L=/tmp/ed_$tag
  timeout 2400 /usr/bin/time -v env OMP_NUM_THREADS=1 SCT_RSWEEP=1 \
    SCT_RMIN=$rmin SCT_RMAX=$rmax SCT_SMAX=$sm \
    $RN $GD/$g.edges $rmin $((rmin+1)) >/dev/null 2>$L; ex=$?
  echo "  [$tag] $g plane r=$rmin..$rmax smax=$sm -> $(tsec $L)s (exit=$ex)"
  grep -E 'PATTERN-EXPLOSION' $L 2>/dev/null | sed 's/^/      /'
}

echo ""; echo "##### A. COMPRESSION CONTRAST (fixed-r sweep -- reaches the counters) #####"
for r in 3 4 5; do
  echo "--- r=$r : com-amazon (WE WIN 5.0x) ---"; sweep amz_r$r com-amazon.ungraph $r 8
  echo "--- r=$r : email-Eu-core (WE LOSE 19x) ---"; sweep eml_r$r email-Eu-core $r 8
done

echo ""; echo "##### B. email plane: which r dominates? (smax=8) #####"
for r in 3 4 5; do plane epl_r$r email-Eu-core $r $r 8; done

echo ""; echo "##### C. email plane: where does s explode? (r=3..5) #####"
for sm in 5 6 7 8; do plane epl_s$sm email-Eu-core 3 5 $sm; done

echo ""; echo "=== EMAIL_DIAG_DONE $(date) ==="
