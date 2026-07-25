> # !! READ THIS BEFORE USING ANY NUMBER BELOW (2026-07-26) !!
> **Every CND-based ratio in the RQ1/RQ2/RQ3 tables of this file predates the serial-CND rule and was
> measured against PARALLEL CND. They are QUARANTINED by SigmodPlus.md §217 and must not go in the
> paper.** §207 proved parallel CND is both slower and fatter on these inputs, so the corrections do
> not even point the same way: serial CND is FASTER on dense-collab cells (ratios get worse for us)
> but COMPLETES where parallel aborted (web-it's "infeasible" becomes a clean, much stronger 742x).
>
> **The current, serial-protocol, same-machine data lives in `../docs/RESULTS.md` and in
> SigmodPlus.md §216 (single-cell table), §218-219 (the full (r,s) grid over the acceptance roster),
> §221-223 (the index).** The direction also changed: the paper is now INDEX-first (§213), so the
> RQ structure below is being superseded, not merely re-measured.
>
> Retained here because the RQ *questions* and the data inventory are still the right skeleton.

# EXPERIMENTS.md -- the paper's experiment section, live status
(Maintained per user directive 2026-07-08. Source of truth for WHAT goes in the paper and WHICH
data already exists. Detailed narrative lives in SigmodPlus.md §139-141; update BOTH on change.)

## Protocol (binding for every number)
> **STILL BINDING, with two additions (2026-07-26).** (a) CND must be run SERIAL
> (`OMP_NUM_THREADS=1`) from the repo root; every pre-2026-07-21 number here violates that.
> (b) Attribute a metric to the PHASE you changed: on MCE-dominated graphs total-time comparison
> manufactures false wins AND false losses from the same noise (§210).

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
> ## !!! EVERY NUMBER IN THIS SECTION IS VOID !!!
> All of its CND columns are **parallel** CND (§217). Serial CND is faster on dense-collab cells
> and COMPLETES where parallel aborted, so the corrections do not even share a direction.
> The `>=300GB abort` markers in particular are wrong: serial CND FINISHES web-it (3,4) at 101GB
> in 93 min, which is a far stronger and defensible **742x / 226x** (§216).
> **REPLACEMENT: SigmodPlus §216 (single cell) and §218-219 (the full (r,s) grid).**
> Do not re-measure this table; the grid already supersedes it on a bigger, cleaner roster.

Status: REBUILT on UNIFORM r=4, s=5..8 (all graphs same setting, §145). SigmodPlus §145.
Paper Exp-1 = uniform table (Table tab:main). Stable advantage = 8 clique-structured graphs
(5 CND-infeasible + raefsky3 1600x/pkustk11 770x/dblp 35x/nasasrb 15x). Honest: astro 1.5x,
yt parity; hepph/epin DROPPED from headline (pattern wall at r=4, per user 'only need 8 graphs').
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
HONEST MEMORY ROWS: on the social graphs our sweep RSS exceeds CND's (epin 26.2GB vs ~2.0GB,
yt 7.8GB vs ~3.9GB; residue machinery constants); the memory win is at the infeasibility
frontier, not universal -- state it plainly in RQ1 prose.
Alternates: sc-msdoor (infeasible), if a reviewer objects to any pick.
Data: tods2:/home/wenqianz/nsi_main_table/ + /home/wenqianz/nsi_e3/ + /home/wenqianz/nsi_roster/.
TODO: copy the three tods2 result dirs into paper_data/ (git, -f) before writing.

## RQ2: the chain certificate's power -- SELF-comparison (sweep marginal vs native per-cell)
> **STILL VALID as data (ours-vs-ours, no CND involved), but DEMOTED in role.** Per the standing
> instruction that the paper's numbers are vs CND, this is an ABLATION, not a result. It also
> predates the §210 optimization stack, so the absolute times are stale even though the ratios
> illustrate the certificate correctly.

Status: DATA COMPLETE (single-trial). SigmodPlus §139 (E1 Phase B).
Marginal-cell headlines (same machine): dblp5 s=9: 0.15s vs >7200s (>48000x); dblp4 s=8: 0.07s
vs 655.5s (9364x); hepph s=7: 10.4s vs >7200s (>692x); astro s=8: 10.7s vs >7200s (>673x);
astro s=7: 3.3s vs 1466s (439x); hepph s=6: 2.2s vs 914s (410x). Spectrum totals: hepph >50.7x,
astro >100x, dblp4 74.5x, dblp5 >258x, webit 5.3x, epin 1.43x, yt 1.61x.
KEEP STRICTLY SEPARATE from RQ1 in the paper (self vs competitor comparison).

## RQ3: the index -- build/size/query
> ## !!! SUPERSEDED -- and the old numbers made the index look WORSE than it is !!!
> The table below measures the FAT index format. Its 40.30 B/r-clique on epin and 34.23 on yt are
> **larger than the archive they replace** (4r + 4*cells = 28 B/r-clique at r=3 over 4 cells).
> That was an ENCODING problem, diagnosed in §221: the pattern table was 86-99% of the file.
> **REPLACEMENT: §223-224. NSI3 drops every certified pattern record and is 8.9x-1359x smaller
> than the full plane index, verified answer-for-answer on 1.2M spectrum rows across 8 graphs.**
> E4 (the query baseline) is no longer 'PENDING with no tool': `region_native/nsi_baseline.cpp`
> builds the sorted archive and probes it. Pilot: index 5.5x smaller AND 1.4x faster than the probe.

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

## RQ-SCAL: scalability with (r,s) -- NEW, the strongest experiment (SigmodPlus §149)
> **STILL VALID (CND-independent) and now EXPLAINED.** The 'advantage grows with r' observation
> has a mechanism and a formula: the class alphabet is r-INDEPENDENT while mult(P) grows
> binomially, so compression = E_P[prod_c C(n_c,b_c)] rises with r (§220), and W = hostSz/compression
> falls monotonically in r on every graph measured (§209b). Cite the mechanism, not just the curve.

Status: DATA COMPLETE (paper_data/scalability/scal_2026-07-09.tsv). Paper Exp-2/Exp-3 + 2 figures.
Panel A (vary r=3..7): webit 0.8/1.5/2.9/5.3/15s, dblp 2.9/5.9/18/47/97s -- advantage grows with r.
Panel B (fix r=4, smax=6..12): webit 1.41->1.73s, dblp 5.27->6.09s, mem FLAT -- spectrum width is
nearly FREE (each cell ~0.05-0.1s). THE headline structural result, CND-independent.
Panel C: CND single cell (3,4) raefsky3 10.7s/pkustk11 67s/nasasrb 28s/dblp 7.7s vs SpecND whole
r=3 spectrum -- crossover. Figures TODO.

## RQ4: certification anatomy (WHY it works) -- figure
> **STILL VALID, and now load-bearing for the INDEX, not just for speed.** The certified fraction
> is exactly the fraction of pattern records NSI3 can drop (measured 99.99-100% on the roster,
> §223-224), so this figure doubles as the index-size explanation.

Status: data embedded in every sweep log ([nsi-cell] lines); AGGREGATION SCRIPT TODO.
Certified% per cell: FEM/CFD 100% (residue 0; pwtk ~4k/3.86M); web-it 100%, web-uk ALL-MERGEABLE
(closed form, no patterns at all); dblp 100.00% (residue literally 0); collab 99.96-99.97%;
yt 79-81%; epin 47-51%. Figure: per-graph-family certified fraction per cell (bar/heatmap) +
the T5 equality rates (GrQc/HepPh 100%, epin 23.4%, §133) as the diagonal analogue.

## RQ5: honest boundary
> **SUPERSEDED BY A STRONGER, PREDICTIVE VERSION.** 'Social graphs only get 1.4-5.9x' is now a
> computable a-priori statement: W = hostSz_avg / compression, with compression = E_P[prod C(n_c,b_c)],
> both available from the front end alone (`SCT_W_ONLY=1`). There are THREE distinct boundary types,
> not one (§209, §218): compression ~ 1 (email/pokec, twin-free), host-multiplicity smear (ca-HepPh),
> and the pattern wall (dblp-coauthor). Write the characterization, not the apology.

Status: DATA COMPLETE. Social graphs (epin/yt) ~50-80% certified -> 1.4-5.9x only; astro/hepph
parity/loss vs CND (s-flat counting); P10 no-free-spectrum (sketch -- MUST either write the
gadget out or demote before submission, §140 P7). Band/diagonal wall (§131-132) = one scope
paragraph, not a contribution.

## Correctness gates (not a paper experiment; reviewer-facing reproducibility)
- Sweep vs native distribution: bit-exact on every feasible cell (GrQc/HepPh/AstroPh/Epinions +
  E1 phase gates). Deep cells (native exceeds budget): gate N/A, proofs carry.
- Index: 4.5M point queries vs REF dumps, exact (3 graphs incl heavy-residue).
- T5: 0 violations across 156M+ 4-cliques.
- CND precondition: binary has the 8ea7546 fix; PIVOTER_COMPARE spot-check PASSED (hepph 3,5:
  "Optimized vs Reference correctness verified (exact)", tods2:/home/wenqianz/nsi_e3/compare_hepph35.log).

## Remaining runs (from §140, updated)
> **THIS LIST IS OBSOLETE.** It was written for the decomposition-speed paper. The current
> ranked gaps are at the bottom of this file under CURRENT STATUS, and the live checklist is
> `TODO.md` (INDEX TRACK section).

- [ ] E2 multi-trial: >=3 runs/median for RQ1+RQ2+RQ3 headline numbers (script the reruns).
- [ ] E4 query baseline: sorted-table probe vs nsi_query, same-machine.
- [ ] E5b: indexes for the 8 roster graphs (sweep with SCT_INDEX_OUT; webit/webuk/FEM row).
- [ ] RQ4 aggregation script (parse sweep logs -> certified% CSV -> figure).
- [ ] Copy tods2 result dirs into paper_data/ (git -f) + verify compare_hepph35.
- [ ] (stretch) E8 one extra large graph (as-skitter / soc-pokec r=3).

## Data inventory
> **The tods2 paths below still exist but hold PARALLEL-CND era output (§217).** Current output:
> `/data/wenqianz/{grid_scout_batch1,grid_scout_fem,cnd_vs_stack,nsi3}.out` and the indexes in
> `/data/wenqianz/nsi3idx/`.

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

---

## CURRENT STATUS (2026-07-26) -- what actually exists now

### Acceptance standard (user, binding)
>= 8 graphs, all million-scale, 3-4 domains, all excellent. **MET** on 8 graphs / 4 domains, plus two
59M-162M-edge FEM additions. Billion-scale web (it-2004 1.03B edges, uk-2005 783M) is downloaded and
converted but blocked at MCE (>3h budget) -- recorded as an honest scaling boundary, not hidden.

### The roster and its (r,s) grid vs SERIAL CND (§218, §219)
| domain | graph | edges | grid verdict |
|---|---|---|---|
| web | web-it-2004 | 7.2M | CND 0/9 cells feasible (budget or 300GB kill); our whole rows 0.8-3.2s |
| web | web-Google | 5.1M | 9/9 wins, 1.43x -> 131x along r and s |
| collab | com-dblp | 1.05M | all wins, to 626x |
| co-purchase | com-amazon | 926k | all wins, 14.9x -> 163x |
| FEM | pwtk | 5.71M | 9/9, deep cells CND-infeasible |
| FEM | pkustk13 | 3.26M | 9/9, to 1448x |
| FEM | pkustk11 | 2.57M | 9/9, to 2781x |
| FEM | nasasrb | 1.31M | 9/9, to 1661x |
| FEM (large) | Flan_1565 | 59M | CND infeasible at every cell; our r=3 row 282s |
| FEM (large) | Queen_4147 | 163M | CND infeasible at every cell; our r=3 row 895s |
Single-cell headline (§216, 3-trial medians, both sides serial): **web-it-2004 (3,4) 5581s / 101GB
for CND vs 7.52s / 458MB for us = 742x time, 226x memory.**

### Deliberate exclusions, each explained by the characterization (§209, §218, §220)
com-youtube (r=3 losses), soc-pokec (compression 1.001x), com-lj / com-orkut (same twin-free
structure), dblp-coauthor (a THIRD boundary type: the pattern wall -- MCE finishes but incidences
blow the 500M cap). Social is absent by structure, and W predicts that before any build.

### The index (§221-223) -- now the paper's main axis
NSI3, the slim plane index, stores only classOf + class-to-region profiles + the exceptions the
theory cannot reconstruct. Measured against the full NSI2 plane index, same graphs, gate = every
answer identical:
| graph | NSI2 | NSI3 | shrink | gate |
|---|---|---|---|---|
| pkustk13 | 1.86 GB | **1.4 MB** | **1359x** | PASS (150,000 spectrum rows) |
| nasasrb | 825 MB | **1.5 MB** | **558x** | PASS |
| pwtk | 404 MB | 3.2 MB | 127x | PASS |
| pkustk11 | 91 MB | 0.83 MB | 109x | PASS |
| com-dblp | 238 MB | 4.2 MB | 57x | PASS |
| com-amazon | 38 MB | 4.3 MB | 8.9x | PASS |
| ca-CondMat | 18.7 MB | 376 KB | 52x | PASS (120,000 rows, r=3,4,5) |
| ca-GrQc | 4.17 MB | 51.7 KB | 85x | PASS |
Query: cold latency improves on all measured operations (1.4-1.9x, the structure misses cache far
less); warm is at parity after the decreasing-size profile order with early exit; index load 30x
faster. Against a materialized archive at 4r + 4*cells bytes per r-clique, NSI3 is 58x-1100x smaller.

### Remaining experiment gaps (ranked)
1. **E4 query baseline vs the archive on the ROSTER** (the tool exists: region_native/nsi_baseline.cpp
   builds the sorted table and probes it; pilot on ca-GrQc gave index 5.5x smaller AND 1.4x faster
   than the archive probe). Needed for the index paper's query axis.
2. Index build on the roster via the OPTIMIZED path: NSI2/NSI3 are written by the PLANE engine, which
   still has none of the §210 optimizations (debt #1). Either port them or build per-r and merge.
3. Multi-trial medians for every number that will be quoted; the grid used the §fast-protocol
   (1 trial where the ratio is lopsided, refine queue for cells in [1/3, 3]).
4. Billion-scale: the MCE wall, and whether a k-core prefilter is legitimate to report.
