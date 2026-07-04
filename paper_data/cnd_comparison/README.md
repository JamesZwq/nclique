# Existing ours-vs-CND comparison data

Same-machine, fixed-CND (commit 8ea7546), tods2, 1h timeout, OMP=1,
/usr/bin/time -v wall + peak RSS. Provenance: SigmodPlus §101 (2026-06-22).

| File | Grid | Source on server |
|---|---|---|
| `cmp_fixed_2026-06-22.csv` | small/mid graphs, r=2..5 | `/data/wenqianz/cmp_fixed.csv` (driver `cmp_fixed.sh`) |
| `cmp_big_2026-06-22.csv`   | big-clique graphs (HepPh/Stanford/it-2004/pokec) | `/data/wenqianz/cmp_big.csv` (driver `cmp_big.sh`) |

Columns: `graph,r,s,method,wall_s,rss_mb,status`. method in {OURS, CND}.
status: OK / TIMEOUT / ERR_rc137 (OOM-kill) / ERR_rc134 (abort) / ERR_rc1 (load error).

## Headline (from this data)
- OURS WINS high-RS: com-dblp 5,6 70.7s/3.6GB vs CND 1018s/93.8GB (14.4x, 26x leaner, CND near-OOM);
  com-dblp 4,5 2.1x; ca-GrQc 5,6 4.7x/12x; web-it-2004 3,4 8.5s/0.45GB vs CND TIMEOUT/83GB (feasibility win).
- OURS LOST dense: ca-AstroPh 4,6 1479s vs CND 24.7s (60x); ca-HepPh 3,4 1052s vs 11.8s; web-Stanford 3,4 880s vs 57s.

## IMPORTANT TIME-VALIDITY CAVEAT
The OURS column here is the a_Y binary BEFORE §117-121 (it is the SIGMOD-freeze engine).
The CND column is still current (CND has not changed). §117-121 sped the a_Y PEEL up to 6.4x on the
loss cells (e.g. ca-AstroPh 4,6 peel 573s->90s), so the OURS wall for build+peel loss cells is now
PESSIMISTIC. To quote a current ratio, RE-RUN OURS ONLY (region_native/region_native_sct_peel at
main) on the loss cells and pair it with the CND column already here. See SigmodPlus §122 §4.
