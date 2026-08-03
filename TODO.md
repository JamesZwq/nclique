# NSI / SIGMOD paper TODO (mirror of SigmodPlus.md §140 -- that section is the source of truth)

## Experiments (tods2, serial, same-machine only)
- [ ] E1 main table 7 configs (RUNNING on tods2; poller will report)
- [ ] E2 multi-trial (>=3 runs, median) for every paper number (--trials flag)
- [ ] E3 CND same-machine spectrum baseline (cell list in §140.1; neutral OOM markers)
- [x] E4 same-machine query bench (NSI vs materialized-archive probe) -- **DONE, §233**
- [ ] E5 server-scale index build/size/latency for all 7 configs
- [ ] E6 certification-anatomy figure (cert% per cell per graph)
- [x] E7 certificate ablation **DONE in substance (§236)**: no-T3 = **22.6x** slower (cells 67-5406x),
      no-T5 = 2.4-3.0x, answers byte-identical on both graphs; pkustk13 no-T3 completing
- [ ] E8 (stretch) one extra large graph

## Cost anatomy / turn the loss cells around (§208, LIVE -- this is the active thread)
- [ ] G1 probe M = classes-per-leaf (incidence-weighted). Gate: median M/r < 3 => kill deconvolution
- [ ] G2 prototype deconvolution supInit on ca-HepPh (3,5) behind a flag; gate = max rel. error vs the
      current double DP over ALL patterns + supInit-segment wall time. NOT bit-identical: build the
      tolerance (or exact-integer) gate FIRST, the existing bit-exact gates do not apply
- [ ] G3 addDelta grouping on **ca-AstroPh (4,6)** (61.3% of total there), NOT on ca-HepPh (19.7%)
- [ ] G4 recompute W on the 7 win-hunt graphs + the loss set; confirm the threshold. This table is
      publishable on its own as the honest "when does this method apply" characterization
- [x] ~~SLP step 2 grouped walk~~ **DE-GREENLIT (§208)**: §203's 71-163x lumped ungroupable
      bounds-rejects in; the clean redundancy is 2.2-9.8x

> Before adding anything to this list, check `DO_NOT_REPEAT.md`: TWELVE directions are already dead
> with the measurement that killed them, and sixteen traps already produced wrong results.

## INDEX TRACK (§213 decided the paper is index-first; this is now the main line)
- [x] §221/§222 byte anatomy + slim design VERIFIED (cP recoverable, 0 mismatches)
- [x] §223 NSI3 shipped: 52-1359x smaller, bit-identical, faster cold, warm at parity
- [x] **§225 E4 query baseline built** -- `nsi_query INDEX archive` / `archive-bench R QUERIES`,
      plane-aware, in the SAME binary as the index so both are timed identically. Every choice is
      made in the ARCHIVE's favour (see §225), so the gap is a lower bound. `scripts/e4_archive_tods2.sh`
- [x] §226 anatomy of what NSI3 still keeps (`nsi_query INDEX anatomy`, prices a packed encoding
      without writing anything) -> §227 NSI4 packed format: **3.85x** on ca-GrQc, load 1.7x faster,
      60,000 rows byte-identical vs BOTH NSI3 and NSI2. `scripts/nsi4_tods2.sh`
- [x] §228 the plane BUILD: certify before building the leaf maps; skip them, and the per-cell leaf
      view, when a column/cell has no residue. Byte-identical index; ca-GrQc memory 1.40x.
      `scripts/plane228_tods2.sh` for the roster A/B
- [x] all three roster runs DONE and written up: **§233** (query axis), **§234** (NSI4),
      **§235** (the build A/B: faster AND leaner on 8/8, index byte-identical on 8/8)
- [ ] **E2 multi-trial medians** -- now the single most load-bearing gap. §234 caught the same cell
      re-measuring 199 -> 305 ns hours apart on an idle machine. Every latency and every small build
      ratio in §233/§234/§235 is ONE run. The byte-identical gates are not affected.
- [ ] a motivating query workload (currently uniform over r-cliques via `sample --by-clique`); the
      application is interactive multi-resolution dense-subgraph exploration
- [ ] ~~port the §210 stack into the plane engine (old debt #1)~~ **RE-RANKED by §228**: the §210
      stack is a PEEL optimization and the plane build's peel is 13-50%, while the leaf maps were
      28-52% and largely wasted. Fix the waste first (done); re-measure before porting anything.
- [ ] (ranked last, §229) T6 host-1 at r=rMin: ~44% of the first row, ~0% above it, ~4% of the build

## ACCEPTANCE STANDARD (user, 2026-07-23 -- binding)
>= 8 graphs, all million-scale+ (aim for one billion-scale), from 3-4 common domains, ALL excellent
under the index-first metrics (build feasibility, index size, query). Consequences:
- [ ] retire sub-million graphs from the headline roster (raefsky3 733k, ca-* series, epin)
- [ ] W-prescreen the big candidates: dblp-coauthor 30.8M, soc-pokec 30M, com-lj 34M, com-orkut 117M (RUNNING)
- [ ] acquire billion-scale WEB graphs (our strongest domain): LAW it-2004 full 1.15B / uk-2005 936M / sk-2005 1.9B
- [ ] acquire larger FEM matrices (SuiteSparse: nlpkkt160/200/240, Queen_4147, Flan_1565)
- [ ] com-friendster 1.8B (server) -- W-prescreen before committing (social, weak domain)

## Paper
- [ ] P1 verify SIGMOD/VLDB round dates (do NOT trust memory)
- [ ] P2 write in order: Intro+Fig1 -> Theory -> Algorithm -> Index -> Experiments -> Related -> Abstract
      **The narrative architecture is DECIDED: §237** (the spine sentence, the question chain, the
      roadmap table, the 4-term budget, the 4-line kernel shown first, Fig 1's three panels).
      Write TO that section; do not re-derive the story.
- [ ] P6 novelty reading: Burkhardt-Faber 1806.05523, Frohmader flag-KK, ICDE'21 2011.00749
- [ ] P7 P10 gadget: write out formally OR demote to remark (no "sketch" in contributions)
- [ ] P8 create paper dir from acmart template (tracked real dir); commit every edit
- [ ] P9 reproducibility appendix from §134

## Style discipline (binding, from memory)
no em-dashes | contributions = final deliverables | no self-exposed weaknesses |
spectrum-vs-spectrum only | storytelling voice | teacher-paper architecture | 1 sentence/line LaTeX

## Idea backlog (do not block the paper)
P10 gadget theorem | general diagonal transfer | band sub-quotient | T5 engine |
dynamic NSI | parallel sweep | index format v2 | T3-tightness characterization
