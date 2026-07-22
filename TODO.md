# NSI / SIGMOD paper TODO (mirror of SigmodPlus.md §140 -- that section is the source of truth)

## Experiments (tods2, serial, same-machine only)
- [ ] E1 main table 7 configs (RUNNING on tods2; poller will report)
- [ ] E2 multi-trial (>=3 runs, median) for every paper number (--trials flag)
- [ ] E3 CND same-machine spectrum baseline (cell list in §140.1; neutral OOM markers)
- [ ] E4 same-machine query bench (NSI vs sorted-table probe) -- local "parity" claim void until then
- [ ] E5 server-scale index build/size/latency for all 7 configs
- [ ] E6 certification-anatomy figure (cert% per cell per graph)
- [ ] E7 NOCERT ablation (certificate contribution vs shared build)
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
