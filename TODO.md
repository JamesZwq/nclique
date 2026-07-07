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
