# S Trees Over Chains: Measured Bytes And Query Latency

Date: 2026-09-19. Same program, harness and protocol as stage 2
([RESULTS_INDEX.md](RESULTS_INDEX.md)); the only change is the class
definition: hierarchy-equivalence chains ([CHAINS.md](CHAINS.md)) in place
of closed-neighbourhood twins. Evidence: `index_chains.json`,
`index-logs_chains/`. The twin-class run of stage 2 is untouched.

## 1. Problem Summary

Measure, not project, the index of [CHAINS.md](CHAINS.md) Section 4: the
S-tree baseline and the skyline index built over chain classes, against
the stage-2 twin-class layouts, on the five local full-range inputs.

## 2. Current Baseline And Bottlenecks

Baseline: stage-2 S trees over twin classes (the fastest exact
representation found so far). Also reported: the naive representation
(per-vertex S trees plus one explicit value per active (vertex, size)).

## 3. Candidate Algorithm Ideas

S trees over chains; skyline over chains; label alignment (bitmap map),
the last one computed, not built.

## 4. Chosen Approach And Rationale

`index.cpp --graph <path> chains`: build all per-size trees first, group
vertices by the tuple of their own nodes, then run the unchanged layout,
selftest and measurement code with chains as classes. Vertices with no
clique form one class. Build-time checks: every member of a class has
the same own node at every size and the same core value at every size.

## 5. Complexity Discussion

Chain grouping is one map insertion per vertex with a key of length
omega(v) - 1. Everything else is as in stage 2. Whole runs took 2 s
(GrQc) to 334 s (amazon); the community batches dominate.

## 6. Implementation Summary

`index.cpp` (class mode argument, both modes in the selftest),
`run_index.py <mode>` (mode-specific evidence files), `chains.cpp` (the
counts of CHAINS.md).

## 7. Correctness Validation Summary

| Check (both class modes, every graph of the stage-2 selftest set) | Release | ASan/UBSan |
|---|---:|---:|
| Graphs | 34,075 | 34,075 |
| Community queries against brute-force nuclei (two designs, two modes) | 1,210,952 | 1,210,952 |
| Membership answers against brute force | 7,303,180 | 7,303,180 |
| Value answers against the core matrix | 2,879,668 | 2,879,668 |
| `count --selftest` unchanged | 34,075 / 844,230 / 93,057 | same |

On the five real graphs every timed query (41,000 community, 20,001
membership per graph) was answered identically by both designs; the
build-time class checks passed; the core matrix equals the frozen control.

## 8. Experimental Setup

As in stage 2 (same seed, same query mix: own / half / root levels,
20,000 / 20,000 / 1,000 community queries, 20,001 membership, 200,000
values; one warm-up and five timed passes, median). Commands:
`python3 research/r1_skyline_index_20260918/run_index.py chains` from
the repository root.

## 9. Experimental Results

Bytes with values (Block D), measured layouts; the last column is the
label-aligned projection (the vertex-to-chain map replaced by a bitmap,
a rank directory and one range per chain):

| Graph | classes twins -> chains | (class, s) pairs twins -> chains | S trees, twins | S trees, chains | ratio | skyline, chains | aligned (projection) | ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| ca-GrQc | 4,105 -> 716 | 12,704 -> 2,079 | 256,238 | 103,220 | 2.48 | 114,382 | 65,128 | 3.93 |
| ca-HepPh | 9,559 -> 1,136 | 114,694 -> 8,468 | 1,532,406 | 395,336 | 3.88 | 398,736 | 306,065 | 5.01 |
| com-dblp | 256,983 -> 13,459 | 928,980 -> 70,208 | 16,880,402 | 4,447,786 | 3.80 | 4,558,750 | 2,024,433 | 8.34 |
| web-Stanford | 270,303 -> 17,963 | 1,429,928 -> 151,364 | 26,086,394 | 5,762,930 | 4.53 | 6,327,402 | 3,632,412 | 7.18 |
| amazon0302 | 253,766 -> 40,867 | 787,990 -> 143,567 | 17,079,748 | 5,744,638 | 2.97 | 6,283,334 | 3,860,362 | 4.42 |

Blocks of S trees over chains (bytes): vertex-to-chain map and
chain-to-vertex CSR / nodes / (chain, s) pairs / values:

| Graph | map | nodes | pairs | values |
|---|---:|---:|---:|---:|
| ca-GrQc | 44,804 | 30,340 | 19,500 | 8,576 |
| ca-HepPh | 100,612 | 155,672 | 72,292 | 66,760 |
| com-dblp | 2,590,480 | 951,412 | 615,504 | 290,390 |
| web-Stanford | 2,327,080 | 1,064,780 | 1,282,768 | 1,088,302 |
| amazon0302 | 2,260,360 | 1,302,680 | 1,312,008 | 869,590 |

Against the naive representation (per-vertex S trees plus one explicit
value per active pair): GrQc 372,032 (3.6x, aligned 5.7x); HepPh
7,122,228 (18x, 23x); dblp 32,146,128 (7.2x, 15.9x); Stanford 27,920,220
(4.8x, 7.7x); amazon 15,550,984 (2.7x, 4.0x).

Latency, ns per query, S trees over twins / S trees over chains /
skyline over chains; community outputs are identical in size:

| Graph | own | half | root | membership | value |
|---|---|---|---|---|---|
| ca-GrQc | 2,903 / 574 / 2,170 | 3,449 / 675 / 2,433 | 4,200 / 734 / 2,534 | 4.3 / 4.9 / 20.0 | 4.3 / 5.0 |
| ca-HepPh | 9,668 / 2,116 / 11,041 | 10,014 / 2,213 / 11,254 | 11,129 / 2,189 / 10,697 | 13.0 / 8.1 / 55.9 | 5.4 / 5.2 |
| com-dblp | 370,861 / 42,952 / 127,411 | 442,375 / 47,615 / 143,467 | 498,255 / 58,554 / 157,674 | 14.1 / 8.6 / 29.4 | 4.6 / 4.6 |
| web-Stanford | 163,655 / 27,140 / 171,661 | 209,661 / 34,461 / 197,738 | 256,140 / 40,847 / 219,552 | 89.3 / 28.6 / 54.1 | 8.2 / 5.9 |
| amazon0302 | 223,981 / 133,725 / 312,831 | 326,479 / 191,094 / 479,997 | 352,531 / 206,231 / 507,648 | 9.7 / 6.0 / 24.6 | 8.5 / 6.7 |

Per output vertex at the own level: twins 1.75 / 1.81 / 2.25 / 1.98 /
2.09 ns, chains 0.35 / 0.40 / 0.26 / 0.33 / 1.25 ns.

## 10. Analysis Of Runtime And Memory

- Memory: S trees over chains are 2.5x to 4.5x smaller than the stage-2
  baseline as measured, and the vertex-to-chain map is now 45 to 58
  percent of the total on the three large graphs; aligned labels would
  remove most of it (3.9x to 8.3x, projection). The (chain, s) pairs, the
  block that S trees pay per vertex, shrank 5.5x to 20x; the values block
  shrank 6x to 8.6x because a chain stores one trajectory.
- Speed: community listing is 5x to 8.6x faster on four graphs and 1.7x
  on amazon at identical output, because a community is now a short list
  of chains whose vertex groups are copied as contiguous blocks (0.26 to
  0.40 ns per output vertex, memory-copy speed) instead of one class per
  vertex; amazon's chains are small (6.4 vertices on average), hence the
  smaller gain. Membership is 1.5x to 3.1x faster on the three large
  graphs (fewer distinct classes, better locality) and within noise on
  the two small ones. Values are unchanged.
- The skyline index over chains is larger than S trees over chains on
  every graph and 2x to 5x slower: with chains, the per-pair block is
  small and the skyline's node, location and cross-size blocks dominate.
  The skyline dedup is unnecessary once chains are used.

## 11. Failure Cases / Limitations

- The aligned-label column is a projection; building it requires the
  index to relabel vertices (a permutation, 4 bytes per vertex, unless
  the graph is stored in that order).
- amazon0302 gains least on both axes (small chains).
- Single machine, warm cache, in-memory; medians of five passes.
- chains <= canonical nodes is observed on all five inputs, not proved.
- Only r = 1.

## 12. Final Conclusion

Chains give what the skyline could not: an index that is both smaller
(2.5x to 4.5x measured, 3.9x to 8.3x projected with aligned labels; 2.7x
to 18x against the naive representation) and faster (community listing
1.7x to 8.6x, membership up to 3.1x) than the best previous exact
representation, with correctness checked against brute force on 34,075
graphs and against the twin-class designs on every timed query.

## 13. Next Recommended Improvements

1. Build the aligned-label layout and measure it (the only unmeasured
   column).
2. Prove or refute chains <= canonical nodes; it bounds the index by the
   hierarchy size rather than by n.
3. Extend chains to r >= 2 (edges and r-cliques as the items), where the
   per-item collapse should be larger still.
