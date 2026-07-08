# EXPERIMENTS.md -- the paper's experiment section, live status
(Maintained per user directive 2026-07-08. Source of truth for WHAT goes in the paper and WHICH
data already exists. Detailed narrative lives in SigmodPlus.md §139-141; update BOTH on change.)

## Protocol (binding for every number)
- SAME-MACHINE ONLY: all paper numbers from tods2 (96-core, 503GB; serial single-thread runs;
  another user's idle IDE daemons noted for the reproducibility paragraph). Local M-series
  numbers are sanity checks only, NEVER in tables.
- Guards: natives/CND 2h budget per cell (timeout -> "exceeds budget" marker), CND 300GB
  prlimit (rc134 -> "memory budget exceeded" marker). Neutral marker semantics in prose.
- Wrapping: /usr/bin/time -v (wall + peak RSS) on every run.
- Correctness: exact algorithm + proofs (docs/nsi_theorems.md); bit-exact gates vs the native
  engine wherever the native cell is feasible; NO "correctness experiments" in the paper.
- PENDING protocol upgrade: E2 multi-trial (>=3 runs, median) before numbers are camera-ready.

## RQ1: spectrum cost vs the per-cell baseline (CND) -- THE MAIN TABLE
Status: DATA COMPLETE (single-trial). SigmodPlus §139 (E1/E3) + §141 (roster).
The 8 big-advantage graphs (picked §141, domain-balanced):
  graph          row      NSI sweep (RSS)   CND spectrum (RSS)      advantage
  web-it-2004    3,4-7    8.1s  (0.5GB)     4x >=300GB abort        infeasible vs 8.1s
  web-uk-2005    3,4-7    7.5s  (0.37GB)    4x >=300GB abort        infeasible vs 7.5s (ALL-CLOSED-FORM)
  sc-pkustk13    4,5-8    81s   (29.4GB)    4x >=300GB abort        infeasible vs 81s
  sc-pwtk        4,5-8    22s   (7.6GB)     4x >=300GB abort        infeasible vs 22s
  sc-ldoor       4,5-8    18s   (7.7GB)     4x >=300GB abort        infeasible vs 18s
  raefsky3       4,5-8    0.6s  (0.2GB)     306.7s (27.9GB)         511x time, 140x mem
  sc-pkustk11    4,5-8    3.7s  (1.4GB)     1349.6s (155.6GB)       365x time, 111x mem
  com-dblp       5,6-9    28.3s (2.7GB)     4261s (257.3GB)         150.6x time, 95x mem
Honest secondary rows (same table or adjacent): dblp r4 22.8x; nasasrb 16x; epin 5.9x; yt 5.8x;
astro 1.65x; hepph 0.82x (CND WINS -- keep, framed by CND's s-flat combinatorial counting).
Alternates: sc-msdoor (infeasible), if a reviewer objects to any pick.
Data: tods2:/home/wenqianz/nsi_main_table/ + /home/wenqianz/nsi_e3/ + /home/wenqianz/nsi_roster/.
TODO: copy the three tods2 result dirs into paper_data/ (git, -f) before writing.

## RQ2: the chain certificate's power -- SELF-comparison (sweep marginal vs native per-cell)
Status: DATA COMPLETE (single-trial). SigmodPlus §139 (E1 Phase B).
Marginal-cell headlines (same machine): dblp5 s=9: 0.15s vs >7200s (>48000x); dblp4 s=8: 0.07s
vs 655.5s (9364x); hepph s=7: 10.4s vs >7200s (>692x); astro s=8: 10.7s vs >7200s (>673x);
astro s=7: 3.3s vs 1466s (439x); hepph s=6: 2.2s vs 914s (410x). Spectrum totals: hepph >50.7x,
astro >100x, dblp4 74.5x, dblp5 >258x, webit 5.3x, epin 1.43x, yt 1.61x.
KEEP STRICTLY SEPARATE from RQ1 in the paper (self vs competitor comparison).

## RQ3: the index -- build/size/query
Status: E5 DATA COMPLETE; E4 (baseline) PENDING. SigmodPlus §139 (E5 table).
  config  index-MB  B/r-clique  load   point-q  spectrum-q     (200k random r-clique queries)
  hepph   33.5      10.47       0.62s  200ns    354ns
  astro   142.4     15.60       2.44s  371ns    574ns
  dblp4   34.5      2.17        0.61s  150ns    275ns
  dblp5   95.7      0.38        1.70s  130ns    253ns
  webit   14.5      0.045       0.26s  81ns     154ns
  epin    59.2      40.30       0.88s  299ns    421ns
  yt      83.7      34.23       1.41s  215ns    328ns
Headlines: webit 0.045 B/r-clique for the WHOLE spectrum (raw per-cell table would be ~27GB,
~1868x); dblp5 0.38B (~307x); index write adds ~0 to the sweep. Correctness: 4.5M point queries
exact vs REF (§136 gates, local -- rerun gates on tods2 optional).
E4 TODO: same-machine sorted-table probe baseline (build lookup tables from REF dumps where
writable, compare per-probe latency + table size + build cost). Roster graphs' indexes TODO
(only the 7 E5 configs have .nsi files: tods2:/home/wenqianz/nsi_e5/).

## RQ4: certification anatomy (WHY it works) -- figure
Status: data embedded in every sweep log ([nsi-cell] lines); AGGREGATION SCRIPT TODO.
Certified% per cell: FEM/CFD 100% (residue 0; pwtk ~4k/3.86M); web-it 100%, web-uk ALL-MERGEABLE
(closed form, no patterns at all); dblp 100.00% (residue literally 0); collab 99.96-99.97%;
yt 79-81%; epin 47-51%. Figure: per-graph-family certified fraction per cell (bar/heatmap) +
the T5 equality rates (GrQc/HepPh 100%, epin 23.4%, §133) as the diagonal analogue.

## RQ5: honest boundary
Status: DATA COMPLETE. Social graphs (epin/yt) ~50-80% certified -> 1.4-5.9x only; astro/hepph
parity/loss vs CND (s-flat counting); P10 no-free-spectrum (sketch -- MUST either write the
gadget out or demote before submission, §140 P7). Band/diagonal wall (§131-132) = one scope
paragraph, not a contribution.

## Correctness gates (not a paper experiment; reviewer-facing reproducibility)
- Sweep vs native distribution: bit-exact on every feasible cell (GrQc/HepPh/AstroPh/Epinions +
  E1 phase gates). Deep cells (native exceeds budget): gate N/A, proofs carry.
- Index: 4.5M point queries vs REF dumps, exact (3 graphs incl heavy-residue).
- T5: 0 violations across 156M+ 4-cliques.
- CND precondition: binary has the 8ea7546 fix; PIVOTER_COMPARE spot-check launched (hepph 3,5,
  tods2:/home/wenqianz/nsi_e3/compare_hepph35.log) -- CHECK RESULT.

## Remaining runs (from §140, updated)
- [ ] E2 multi-trial: >=3 runs/median for RQ1+RQ2+RQ3 headline numbers (script the reruns).
- [ ] E4 query baseline: sorted-table probe vs nsi_query, same-machine.
- [ ] E5b: indexes for the 8 roster graphs (sweep with SCT_INDEX_OUT; webit/webuk/FEM row).
- [ ] RQ4 aggregation script (parse sweep logs -> certified% CSV -> figure).
- [ ] Copy tods2 result dirs into paper_data/ (git -f) + verify compare_hepph35.
- [ ] (stretch) E8 one extra large graph (as-skitter / soc-pokec r=3).

## Data inventory
  paper_data/cnd_comparison/           §124 old per-cell grids (context only; superseded for
                                       spectrum claims by E3)
  paper_data/diag_{astro,dblp}_baseline_2026-07-07.tsv   §131 diagonal U-shape
  paper_data/band_{astro,dblp}_2026-07-07.tsv            §132 band engine
  tods2:/home/wenqianz/nsi_main_table/ E1 (sweeps + natives, v2_summary.log)
  tods2:/home/wenqianz/nsi_e3/         E3 CND cells (e3_summary.log, cnd_*.log, compare log)
  tods2:/home/wenqianz/nsi_e5/         E5 indexes (.nsi) + query workloads + e5_summary.log
  tods2:/home/wenqianz/nsi_roster/     §141 roster (sweep_*.log, cnd_*.log, roster.log)
  docs/nsi_theorems.md                 the theory section source
  SigmodPlus.md §128-141               full narrative + every table above
